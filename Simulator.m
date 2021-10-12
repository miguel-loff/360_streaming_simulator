function [final_qoe, final_qoe_adjusted] = Simulator(nr_users, buffer_request_threshold, scenario, latency, debug_mode, video_allocation, trajectory_allocation, cqi_allocation, initial_buffering_time)
    %%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%
    
    % Simulation settings
    sim.debug_mode = debug_mode; % if true, it will calculate the cumulative QoE every ms and plot the graphs
    sim.t_max = 3 * 60 * 1000; % simultation time - in ms
    sim.nr_users = nr_users; % number of users
    sim.latency = latency; % request-response latency - in ms (fixed value)
    sim.scenario = scenario; % simulation scenario - monolithic, monolithic_cubemap, viewport_dependent, tiles, tiles_partial, viewport_only
    sim.nr_videos = 3; % number of different videos
    sim.rate_adaptation_mode = 'QAAD'; % rate adaptation algorithm - default, QAAD or QDASH
    sim.scheduling_mode = 'PF2';
    sim.nr_rb_per_tti = 106; % number of resource blocks per TTI - TS 38.101-1 or 38.104, table 5.3.2-1 (106 for SCS 15 KHz, BW 20 MHz)
    sim.nr_bits_per_rb_per_cqi = [44, 68, 109, 174, 253, 340, 427, 553, 695, 789, 960, 1128, 1307, 1478, 1605]; % number of bits per resource block - based on the formula in TS 38.306
    sim.buffer_request_threshold = buffer_request_threshold; % buffer level to trigger the request of a new segment - in ms
    switch sim.scenario % these settings will be different depending on the chosen scenario
        case {'monolithic', 'monolithic_cubemap'}
            sim.minimum_period_to_request = 5; % minimum elapsed time to request a new segment - in ms
            sim.nr_segments_requested_initially = 5; % number of segments to reqeust in the initial buffering
            sim.nr_segments_requested_buffering = 5; % number of segments to request after a stall
            sim.qaad.marginal_buffer_length = floor((0.8 * sim.buffer_request_threshold))/1000; %5 marginal buffer length to improve segment quality in QAAD - in s
            sim.qaad.minimal_buffer_length = floor((0.2 * sim.buffer_request_threshold))/1000; %3 critical buffer length in QAAD - in s
        case {'viewport_dependent', 'tiles', 'tiles_partial'}
            sim.minimum_period_to_request = 5; % minimum elapsed time to request a new segment - in ms
            sim.nr_segments_requested_initially = 1; % number of segments to reqeust in the initial buffering
            sim.nr_segments_requested_buffering = 1; % number of segments to request after a stall
            sim.qaad.marginal_buffer_length = floor((0.8 * sim.buffer_request_threshold))/1000; %0.09 marginal buffer length to improve segment quality in QAAD - in s
            sim.qaad.minimal_buffer_length = floor((0.2 * sim.buffer_request_threshold))/1000; %0.09 critical buffer length in QAAD - in s
        case {'viewport_only', 'viewport_only_modified'}
            sim.minimum_period_to_request = 0; % minimum elapsed time to request a new segment - in ms
            sim.nr_segments_requested_initially = 1; % number of segments to reqeust in the initial buffering
            sim.nr_segments_requested_buffering = 1; % number of segments to request after a stall
            sim.qaad.marginal_buffer_length = floor((0.8 * sim.buffer_request_threshold))/1000; %0.003 marginal buffer length to improve segment quality in QAAD - in s
            sim.qaad.minimal_buffer_length = floor((0.2 * sim.buffer_request_threshold))/1000; %0.003 critical buffer length in QAAD - in s
        otherwise
            error('Scenario is not valid');
    end

    % Video settings
    switch sim.scenario % these settings will be different depending on the chosen scenario
        case {'monolithic', 'monolithic_cubemap', 'viewport_dependent', 'tiles', 'tiles_partial'}
            video.segment_length = 1000; % length for each segment - in ms
        case {'viewport_only', 'viewport_only_modified'}
            video.segment_length = 40; % length for each segment - in ms
        otherwise
            error('Scenario is not valid');
    end
    video.psnr_levels = NaN; % V-PSNR levels corresponding to the  6 quality levels - in dB
    video.representations_bitrate = NaN; % bitrates of the 6 quality levels - in Mbps
    video.qualities.levels = NaN; % interpolated quality levels
    video.qualities.psnr = NaN; % equivalent interpolated V-PSNR - in Mbps
    video.nr_bits_per_segment = NaN; % number of bits for each segment
    video.angle = NaN; % angles corresponding to the angular difference between requested direction and viewing direction
    video.delta_psnr = NaN; % impact in the V-PSNR corresponding to the angular differences
    video.default_orientations = NaN; % store the pre defined possible viewport centers for viewport_dependent, tiles and tiles_partial scenarios
    video = repmat(video,sim.nr_videos,1); % create array of video based on number of videos
    
    % Head movement settings
    trajectory.time = NaN;
    trajectory.angle = NaN;
    trajectory.angle_id = 1;
    trajectory = repmat(trajectory,sim.nr_users,1); % create array of trajectory based on number of users
    
    %%%%%%%%%%%%%%%%%%%%%%% VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
    
    % User variables
    user.cqi_events.time = NaN;
    user.cqi_events.cqi = NaN;
    user.cqi_events.event_id = 1;
    user.buffer = 0; % in ms
    user.cqi = 0;
    user.avg_throughput = 0; % in bps
    user.throughput = 0; % in bps
    user.nr_allocated_rbs = 0;
    user.requested_quality = 1;
    user.requested_direction = [];
    user.nr_received_bits = 0;
    user.nr_target_bits = 0;
    user.nr_complete_segments_in_buffer = 0;
    user.nr_played_segments = 0;
    user.last_request = 0;
    user.requesting_data = false;
    user.requesting_initial_segments = true;
    user.stalled = false;
    user.rebuffering = false;
    user.empty_viewport = false;
    user.start_play_time = 0;
    user.initial_start_play_time = 0;
    user.duration_rebuffering_event = 0;
    user.duration_stalled_event = 0;
    user.duration_initial_buffering = 0;
    user.duration_empty_viewport_event = 0;
    user.nr_rebuffering_events = 0;
    user.nr_empty_viewport_events = 0;
    user.served_qualities = [];
    user.qoe = 0;
    user.adjusted_qoe = 0;
    user.estimated_bandwidth = 0;
    user.video_id = 0;
    user.initial_buffering_time = 1;
    user.angle = 0;
    user = repmat(user,sim.nr_users,1); % create array of user based on number of users

    % create array to store adjusted quality
    adjusted_qualities = NaN(sim.nr_users, sim.t_max);

    % Arrays to store simulation results
    aux_results.buffer = NaN(1, sim.t_max);
    aux_results.cqi = NaN(1, sim.t_max);
    aux_results.throughput = NaN(1, sim.t_max);
    aux_results.qoe = NaN(1, sim.t_max);
    aux_results.estimated_bandwidth = NaN(1, sim.t_max);
    aux_results.requested_quality = NaN(1, sim.t_max);
    aux_results.requested_bitrate = NaN(1, sim.t_max);
    aux_results.adjusted_qoe = NaN(1, sim.t_max);
    aux_results.angular_difference = NaN(1, sim.t_max);
    aux_results = repmat(aux_results,sim.nr_users,1); % create array of results based on number of users
    
    %%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%
    
    % Array to store QoE at the end of the simulation
    final_qoe = NaN(1, sim.nr_users);
    final_qoe_adjusted = NaN(1, sim.nr_users);
    
    %%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%
    
    [video, user] = importVideos(video, sim, user, video_allocation);

    [trajectory, user] = importTrajectories(trajectory, sim, user, trajectory_allocation);

    user = importCQIEvents(user, cqi_allocation);
    
    user = setInitialBufferingTime(user, sim, initial_buffering_time);

    for sim_time = 1:sim.t_max
        user = initialBuffering(user, sim, video, sim_time);
        
        [trajectory, user] = updateTrajectory(user, trajectory, sim_time);

        user = fillBuffer(user, sim, video); % buffer is filled with RBs allocated in previous TTI

        user = calculateEstimatedBandwidths(sim_time, user);
        
        user = updateUserState(user, sim, video, sim_time);
        
        user = drainBuffer(user, video, sim_time);

        user = calculateThroughput(user, sim);
        user = calculateAvgThroughput(user, sim_time);

        [user, adjusted_qualities(:,sim_time)] = calculateAdjustedQuality(user, video, sim, adjusted_qualities(:,sim_time));
        
        if sim.debug_mode == true
            user = calculateQoE(user, video, sim, sim_time, adjusted_qualities);
        end

        for user_id = 1:length(user)
            if user(user_id).nr_complete_segments_in_buffer == 0 && user(user_id).stalled == false && user(user_id).rebuffering == false && user(user_id).requesting_initial_segments == false
                user(user_id).stalled = true;
            elseif user(user_id).stalled == true && user(user_id).requesting_data == false && user(user_id).rebuffering == false && user(user_id).requesting_initial_segments == false
                user(user_id) = rebuffering(user(user_id), sim, video, sim_time);
            elseif user(user_id).requesting_data == false && user(user_id).stalled == false && user(user_id).rebuffering == false && user(user_id).requesting_initial_segments == false
                if user(user_id).buffer <= sim.buffer_request_threshold
                    user(user_id) = rateAdaptation(sim_time, user(user_id), sim, video);
                end
            end
        end

        user = updateCQIs(sim_time, user);

        user = allocateRBs(sim_time, sim, user);
        
        if sim.debug_mode == true
            for user_id = 1:length(user)
                aux_results(user_id).buffer(sim_time) = user(user_id).buffer;
                aux_results(user_id).cqi(sim_time) = user(user_id).cqi;
                aux_results(user_id).throughput(sim_time) = user(user_id).throughput;
                aux_results(user_id).qoe(sim_time) = user(user_id).qoe;
                aux_results(user_id).estimated_bandwidth(sim_time) = user(user_id).estimated_bandwidth;
                aux_results(user_id).requested_quality(sim_time) = user(user_id).requested_quality;
                aux_results(user_id).requested_bitrate(sim_time) = video(user(user_id).video_id).representations_bitrate(user(user_id).requested_quality)*10^6;
                aux_results(user_id).avg_throughput(sim_time) = user(user_id).avg_throughput;
                aux_results(user_id).adjusted_qoe(sim_time) = user(user_id).adjusted_qoe;
                aux_results(user_id).angle(sim_time) = user(user_id).angle;
            end
            
            if mod(sim_time,sim.t_max/100) == 0
                fprintf('Simulation progress: %.0f%%\n',(sim_time/sim.t_max)*100);
            end
        end
    end
    
    if sim.debug_mode == true
        plotGraphs(aux_results);
    end

    user = calculateQoE(user, video, sim, sim_time, adjusted_qualities);
        
    final_qoe = writeFinalQoE(user, final_qoe);
    final_qoe_adjusted = writeFinalQoEAdjusted(user, final_qoe_adjusted);
end

%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%

% Import video values and assign videos to users
function [video, user] = importVideos(video, sim, user, video_allocation)
    for video_id = 1:sim.nr_videos
        switch video_id
            case 1
                file_path_video = strcat(pwd, '/videos/chairliftride.csv');
                file_path_delta_psnr = strcat(pwd, '/delta-psnr/chairliftride.csv');
            case 2
                file_path_video = strcat(pwd, '/videos/skateboardinlot.csv');
                file_path_delta_psnr = strcat(pwd, '/delta-psnr/skateboardinlot.csv');
            case 3
                file_path_video = strcat(pwd, '/videos/kiteflite.csv');
                file_path_delta_psnr = strcat(pwd, '/delta-psnr/kiteflite.csv');
        end

        % import video csv
        data = readtable(file_path_video,'ReadRowNames',true);

        % PSNR levels corresponding to the 7 qualities
        video(video_id).psnr_levels = data.V_PSNR(1:7)';

        % bitrates of the 7 qualities
        video(video_id).representations_bitrate = data.(sim.scenario)(1:7)';
        
        % interpolated quality levels (1:0.1:7)
        video(video_id).qualities.levels = 1:0.1:length(video(video_id).psnr_levels);
        % equivalent V-PSNR corresponding to incremental quality levels
        video(video_id).qualities.psnr = interp1(1:length(video(video_id).psnr_levels), video(video_id).psnr_levels, video(video_id).qualities.levels, 'linear','extrap');
        
        % number of bits for each segment
        video(video_id).nr_bits_per_segment = video(video_id).representations_bitrate * 10^6 * 10^-3 * video(video_id).segment_length;
        
        % import delta psnr csv
        data = readtable(file_path_delta_psnr);
        
        % angles and corresponding delta v-psnr
        video(video_id).angle = data.Angle';
        video(video_id).delta_psnr = data.(sim.scenario)';
        
        % interpolate values for 1 degree intervals
        video(video_id).delta_psnr = interp1(video(video_id).angle, video(video_id).delta_psnr, -180:180, 'spline');
        video(video_id).angle = -180:180;
        
        % normalize delta psnr values
        video(video_id).delta_psnr = video(video_id).delta_psnr - max(video(video_id).delta_psnr);
        
        switch sim.scenario
            case 'viewport_dependent'
                video(video_id).default_orientations = [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180];
            case {'tiles', 'tiles_partial'}
                % these orientations are admitting 4 tiles are always downloaded
                video(video_id).default_orientations = [-180, -135, -90, -45, 0, 45, 90, 135, 180];
            otherwise
                video(video_id).default_orientations = -180:180;
        end        
    end
    
    % randomly attribute one of the videos to each user
    for user_id = 1:sim.nr_users
        user(user_id).video_id = video_allocation(user_id);
    end
end

% Import trajectories
function [trajectory, user] = importTrajectories(trajectory, sim, user, trajectory_allocation)
    for user_id = 1:sim.nr_users
        switch user(user_id).video_id
            case 1
                file_path = strcat(pwd, '/trajectories/chairliftride_user_', int2str(trajectory_allocation(user_id)), '.csv');
            case 2
                file_path = strcat(pwd, '/trajectories/skateboardinlot_user_', int2str(trajectory_allocation(user_id)), '.csv');
            case 3
                file_path = strcat(pwd, '/trajectories/kiteflite_user_', int2str(trajectory_allocation(user_id)), '.csv');
        end
                
        data = readtable(file_path);
        
        trajectory(user_id).time = data.Time';
        trajectory(user_id).angle = data.Angle';
    end
    
    % load first angle (sim_time = 1)
    for user_id = 1:sim.nr_users
        user(user_id).angle = trajectory(user_id).angle(trajectory(user_id).angle_id);
        trajectory(user_id).angle_id = trajectory(user_id).angle_id + 1;
    end
end

% Import CQI events from random file
function user = importCQIEvents(user, cqi_allocation)
    for user_id = 1:length(user)
        file_path = strcat(pwd, '/cqis-events', '/cqi_event-', int2str(cqi_allocation(user_id)), '.txt');
        
        aux = importdata(file_path);
        
        user(user_id).cqi_events.time = aux(:,1)';
        user(user_id).cqi_events.time = 1000 * user(user_id).cqi_events.time; % convert time to ms
        user(user_id).cqi_events.cqi = aux(:,2)';
    end
    
    % load first cqi (time = 0)
    for user_id = 1:length(user)
        user(user_id).cqi = user(user_id).cqi_events.cqi(user(user_id).cqi_events.event_id);
        user(user_id).cqi_events.event_id = user(user_id).cqi_events.event_id + 1;
    end 
end

% Set the time each user begins its initial buffering
function user = setInitialBufferingTime(user, sim, initial_buffering_time)
    for user_id = 1:sim.nr_users
        user(user_id).initial_buffering_time = initial_buffering_time(user_id);
    end
end

% Enter initial buffering state
function user = initialBuffering(user, sim, video, sim_time)  
    if sim_time <= 200
        for user_id = 1:length(user)
            if user(user_id).initial_buffering_time == sim_time
                user(user_id).requested_quality = 1;
                user(user_id).requesting_data = true;
                user(user_id).requesting_initial_segments = true;
                user(user_id).last_request = sim_time;
                user(user_id).nr_target_bits = sim.nr_segments_requested_initially * video(user(user_id).video_id).nr_bits_per_segment(user(user_id).requested_quality);
                
                for segment_id = 1:sim.nr_segments_requested_initially
                    user(user_id).served_qualities = [user(user_id).served_qualities 1];
                    
                    % calculate requested angle
                    [~, best_angle_id] = min(abs(video(user(user_id).video_id).default_orientations - user(user_id).angle));
                    user(user_id).requested_direction = [user(user_id).requested_direction video(user(user_id).video_id).default_orientations(best_angle_id)];
                end
            end
        end
    end
end

% Update the angle value of the trajectory
function [trajectory, user] = updateTrajectory(user, trajectory, sim_time)
    for user_id = 1:length(user)
        if trajectory(user_id).angle_id > length(trajectory(user_id).time)
            return
        end
        
        if trajectory(user_id).time(trajectory(user_id).angle_id) == sim_time
            % if in one of these states, the user doesn't move
            if user(user_id).requesting_initial_segments == false && user(user_id).stalled == false && user(user_id).rebuffering == false
                user(user_id).angle = trajectory(user_id).angle(trajectory(user_id).angle_id);
                trajectory(user_id).angle_id = trajectory(user_id).angle_id + 1;
            else
                % delay the trajectory times by 1
                trajectory(user_id).time(trajectory(user_id).angle_id:end) = trajectory(user_id).time(trajectory(user_id).angle_id:end) + 1;
            end
        end
    end
end

% Fill buffer according to allocated RBs last TTI
function user = fillBuffer(user, sim, video)
    %%% only load the milliseconds into the buffer when the segment is complete
    for user_id = 1:length(user)
        user(user_id).nr_received_bits = user(user_id).nr_received_bits + user(user_id).nr_allocated_rbs * sim.nr_bits_per_rb_per_cqi(user(user_id).cqi);
        
        if user(user_id).nr_received_bits >= user(user_id).nr_target_bits && user(user_id).nr_target_bits ~= 0 && user(user_id).requesting_data == true
            user(user_id).buffer = user(user_id).buffer + video(user(user_id).video_id).segment_length;
        end
    end
end

% Estimate current bandwidth based on sampled throughput
function user = calculateEstimatedBandwidths(sim_time, user)
    weight = 0.3; % between 0 and 1 (0.5 seems to be the value that most closely resembles the original behavior)

    for user_id = 1:length(user)
        if user(user_id).nr_received_bits >= user(user_id).nr_target_bits && user(user_id).nr_target_bits ~= 0
            time_to_download = sim_time - user(user_id).last_request;
            bw_sample = user(user_id).nr_target_bits/(time_to_download/1000);

            if user(user_id).estimated_bandwidth == 0
                user(user_id).estimated_bandwidth = bw_sample;
            end

            user(user_id).estimated_bandwidth = weight * user(user_id).estimated_bandwidth + (1 - weight) * bw_sample;
        end
    end

% Original bandwidth estimation algorithm
% It doesn't work for adaptive scenarios because of the lower segment request frequency (designed for continuous requests)
% If using, add this to the beginning of the code:
% sim.qaad.theta = 0.3; % sampling period for the bandwidth estimation algorithm - in s
% sim.qaad.w = 0.875; % weight factor for the moving average for the bandwidth estimation algorithm
% And replace the function declaration by:
% user = calculateEstimatedBandwidths(sim_time, sim, user, [aux_results(1).throughput; aux_results(2).throughput, ...]);
%
%     for user_id = 1:length(user)
%         % the bandwidth starts being estimated after the initial buffering
%         if user(user_id).requesting_initial_segments == false && sim_time ~= 1
%             % make new estimation every theta seconds
%             if mod(sim_time - user(user_id).duration_initial_buffering - 1, sim.qaad.theta*1000) == 0
%                 bw_sample = 0;
%                 
%                 % if there is still less than theta seconds in the simulation, sample all the way back
%                 if sim_time <= sim.qaad.theta * 1000
%                     for time_id = 1:sim_time-1
%                         bw_sample = bw_sample + aux_results(user_id,time_id)/1000;
%                     end
%                     bw_sample = bw_sample/((sim_time-1)/1000);
%                 else
%                     for time_id = sim_time-(sim.qaad.theta*1000):sim_time-1
%                         bw_sample = bw_sample + aux_results(user_id,time_id)/1000;
%                     end
%                     bw_sample = bw_sample/sim.qaad.theta;
%                 end
%                 
%                 if user(user_id).estimated_bandwidth == 0
%                     user(user_id).estimated_bandwidth = bw_sample;
%                 end
%                 
%                 user(user_id).estimated_bandwidth = sim.qaad.w * user(user_id).estimated_bandwidth + (1 - sim.qaad.w) * bw_sample;
%             end
%         end
%     end
end

% Drain buffer
function user = drainBuffer(user, video, sim_time)
    for user_id = 1:length(user)
        if user(user_id).stalled == false && user(user_id).rebuffering == false && user(user_id).requesting_initial_segments == false
            user(user_id).buffer = user(user_id).buffer - 1;
            
            %%% meti +1 para funcionar com drainBuffer depois de updateUserState
            if mod(sim_time - user(user_id).start_play_time + 1, video(user(user_id).video_id).segment_length) == 0
                user(user_id).nr_complete_segments_in_buffer = user(user_id).nr_complete_segments_in_buffer - 1;
                user(user_id).nr_played_segments = user(user_id).nr_played_segments + 1;
            end
        end
    end
end

% Update the user state
function user = updateUserState(user, sim, video, sim_time)
    for user_id = 1:length(user)
        if user(user_id).requesting_initial_segments == true
            user(user_id).duration_initial_buffering = user(user_id).duration_initial_buffering + 1;
            
            if user(user_id).nr_received_bits >= user(user_id).nr_target_bits && user(user_id).nr_target_bits ~= 0
                user(user_id).requesting_initial_segments = false;
                user(user_id).requesting_data = false;
                user(user_id).nr_received_bits = 0;
                user(user_id).nr_complete_segments_in_buffer = user(user_id).nr_complete_segments_in_buffer + sim.nr_segments_requested_initially;
                user(user_id).initial_start_play_time = sim_time;
                user(user_id).start_play_time = sim_time;
            end
        elseif user(user_id).stalled == true && user(user_id).rebuffering == false
            user(user_id).duration_stalled_event = user(user_id).duration_stalled_event + 1;
            
            if user(user_id).nr_received_bits >= user(user_id).nr_target_bits
                user(user_id).requesting_data = false;
                user(user_id).nr_received_bits = 0;
                user(user_id).nr_complete_segments_in_buffer = user(user_id).nr_complete_segments_in_buffer + 1;
                
                if sim.nr_segments_requested_buffering == 0 % user will only wait for the segment to be fully loaded and it will continue without entering rebuffering
                    user(user_id).start_play_time = sim_time;
                    user(user_id).stalled = false;
                    user(user_id).nr_rebuffering_events = user(user_id).nr_rebuffering_events + 1;
                end
            end
        elseif user(user_id).rebuffering == true
            user(user_id).duration_rebuffering_event = user(user_id).duration_rebuffering_event + 1;
            
            if user(user_id).nr_received_bits >= user(user_id).nr_target_bits
                user(user_id).rebuffering = false;
                user(user_id).requesting_data = false;
                user(user_id).nr_received_bits = 0;
                user(user_id).nr_complete_segments_in_buffer = user(user_id).nr_complete_segments_in_buffer + sim.nr_segments_requested_buffering;
                user(user_id).start_play_time = sim_time;
                user(user_id).stalled = false;
            end
        elseif user(user_id).nr_received_bits >= video(user(user_id).video_id).nr_bits_per_segment(user(user_id).requested_quality)
            user(user_id).requesting_data = false;
            user(user_id).nr_received_bits = 0;
            user(user_id).nr_complete_segments_in_buffer = user(user_id).nr_complete_segments_in_buffer + 1;
        end
    end
end

% Calculate bitrate this ms based on allocated RBs
function user = calculateThroughput(user, sim)
    for user_id = 1:length(user)
        user(user_id).throughput = user(user_id).nr_allocated_rbs * 1000 * sim.nr_bits_per_rb_per_cqi(user(user_id).cqi);
    end
end

% Calculate new average throughput based on the new measure
function user = calculateAvgThroughput(user, sim_time)
    for user_id = 1:length(user)
        user(user_id).avg_throughput = (user(user_id).avg_throughput*(sim_time-1) + user(user_id).throughput)/sim_time;
    end
end

% Calculate actual quality based on angular difference between requested direction and viewing direction
function [user, adjusted_qualities] = calculateAdjustedQuality(user, video, sim, adjusted_qualities)
    for user_id = 1:length(user)
        % calculate de angular difference between the direction the user is looking now and the requested direction
        try
            angular_difference = user(user_id).angle - user(user_id).requested_direction(user(user_id).nr_played_segments + 1);
        catch
            angular_difference = NaN;
            continue
        end
        
        while angular_difference > 180
            angular_difference = angular_difference - 360;
        end

        while angular_difference < -180
            angular_difference = angular_difference + 360;
        end

        % calculate the quality impact in dB of this angular difference
        quality_impact = video(user(user_id).video_id).delta_psnr(181 + angular_difference);

        % calculate the actual PSNR
        original_psnr = video(user(user_id).video_id).psnr_levels(user(user_id).served_qualities(user(user_id).nr_played_segments + 1));
        adjusted_psnr = original_psnr + quality_impact;

        % find closest quality level
        [~, adjusted_quality_idx] = min(abs(video(user(user_id).video_id).qualities.psnr - adjusted_psnr));
        adjusted_quality = video(user(user_id).video_id).qualities.levels(adjusted_quality_idx);

        adjusted_qualities(user_id) = adjusted_quality;
    
        % if viewport is not sufficiently covered by an image in tiles_partial or viewport_only, it counts as a stall
        if ((angular_difference < -50 || angular_difference > 50) && sim.scenario == "tiles_partial") || ((angular_difference < -10 || angular_difference > 10) && sim.scenario == "viewport_only") || ((angular_difference < -20 || angular_difference > 20) && sim.scenario == "viewport_only_modified")
            if user(user_id).empty_viewport == false
                user(user_id).nr_empty_viewport_events = user(user_id).nr_empty_viewport_events + 1;
                user(user_id).empty_viewport = true;
            end
            
            user(user_id).duration_empty_viewport_event = user(user_id).duration_empty_viewport_event + 1;
        elseif user(user_id).empty_viewport == true
            user(user_id).empty_viewport = false;
        end
    end
end

% Calculate cumulative QoE up until current ms
function user = calculateQoE(user, video, sim, sim_time, adjusted_qualities)
    for user_id = 1:length(user)
        %%% calculate QoE WITHOUT the effect of the trajectory %%%
        % aggregate all values that influence rebuffering duration
        duration_rebuffering_events = user(user_id).duration_rebuffering_event + user(user_id).duration_initial_buffering + user(user_id).duration_stalled_event;
        
        % calculate mean duration of rebuffering events    
        if user(user_id).nr_rebuffering_events ~= 0
            avg_duration_rebuffering_events = (duration_rebuffering_events/1000)/user(user_id).nr_rebuffering_events;
        else
            avg_duration_rebuffering_events = duration_rebuffering_events/1000;
        end

        %calculate frequency of rebuffering events
        freq_rebuffering_events = user(user_id).nr_rebuffering_events/(sim_time/1000);

        %calculate F_i
        f_i = (7/8) * max([log(freq_rebuffering_events)/6+1 0]) + (1/8) * (min([avg_duration_rebuffering_events 15])/15);

        %calculate average quality
        avg_quality = mean(user(user_id).served_qualities((sim.nr_segments_requested_initially + 1):user(user_id).nr_played_segments)); % only take into account the segments that were actually played
        
        %calculate std. dev. of quality
        std_dev_quality = std(user(user_id).served_qualities((sim.nr_segments_requested_initially + 1):user(user_id).nr_played_segments), 1);

        %calculate qoe
        user(user_id).qoe = 5.67 * (avg_quality/length(video(user(user_id).video_id).representations_bitrate)) - 6.72 * (std_dev_quality/length(video(user(user_id).video_id).representations_bitrate)) + 0.17 - 4.95 * f_i;
        
        if user(user_id).qoe > 5
            user(user_id).qoe = 5;
        elseif user(user_id).qoe < 0
            user(user_id).qoe = 0;
        end
        
        %%% calculate QoE WITH the effect of the trajectory %%%
        if user(user_id).nr_empty_viewport_events > 0
            % aggregate all values that influence rebuffering duration
            duration_rebuffering_events = duration_rebuffering_events + user(user_id).duration_empty_viewport_event;
            nr_rebuffering_events = user(user_id).nr_rebuffering_events + user(user_id).nr_empty_viewport_events;
            
            % calculate mean duration of rebuffering events
            if nr_rebuffering_events ~= 0
                avg_duration_rebuffering_events = (duration_rebuffering_events/1000)/nr_rebuffering_events;
            else
                avg_duration_rebuffering_events = duration_rebuffering_events/1000;
            end
            
             %calculate frequency of rebuffering events
            freq_rebuffering_events = nr_rebuffering_events/(sim_time/1000);

            %calculate F_i
            f_i = (7/8) * max([log(freq_rebuffering_events)/6+1 0]) + (1/8) * (min([avg_duration_rebuffering_events 15])/15);
        end
        
        %calculate new average quality
        initial_time = user(user_id).initial_start_play_time + sim.nr_segments_requested_initially * video(user(user_id).video_id).segment_length + 1;
        avg_quality = mean(adjusted_qualities(user_id, initial_time:sim_time),'omitnan');

        %calculate new std. dev. of quality
        std_dev_quality = std(adjusted_qualities(user_id, initial_time:sim_time), 1, 'omitnan');
        
        %calculate qoe
        user(user_id).adjusted_qoe = 5.67 * (avg_quality/length(video(user(user_id).video_id).representations_bitrate)) - 6.72 * (std_dev_quality/length(video(user(user_id).video_id).representations_bitrate)) + 0.17 - 4.95 * f_i;
        
        if user(user_id).adjusted_qoe > 5
            user(user_id).adjusted_qoe = 5;
        elseif user(user_id).adjusted_qoe < 0
            user(user_id).adjusted_qoe = 0;
        end
    end
end

% Enter rebuffering state
function user = rebuffering(user, sim, video, sim_time)
    user.rebuffering = true;
    user.requested_quality = 1;
    user.requesting_data = true;
    user.nr_rebuffering_events = user.nr_rebuffering_events + 1;
    user.last_request = sim_time;
    user.nr_target_bits = sim.nr_segments_requested_buffering * video(user.video_id).nr_bits_per_segment(user.requested_quality);

    for segment_id = 1:sim.nr_segments_requested_buffering
        user.served_qualities = [user.served_qualities 1];
        
        % calculate requested angle
        [~, best_angle_id] = min(abs(video(user.video_id).default_orientations - user.angle));
        user.requested_direction = [user.requested_direction video(user.video_id).default_orientations(best_angle_id)];
    end
end

% Perform rate adaptation
function user = rateAdaptation(sim_time, user, sim, video)
    switch sim.rate_adaptation_mode
        case 'default'
            if sim_time - user.last_request >= sim.minimum_period_to_request
                if user.avg_throughput > video(user.video_id).representations_bitrate(user.requested_quality) * 10^6 % avg_throughput? 
                    user = requestBetterVideoQuality(user);
                else
                    user = requestWorseVideoQuality(user);
                end

                user.last_request = sim_time;
                user.requesting_data = true;
                user.served_qualities = [user.served_qualities user.requested_quality];
            end
        case 'QAAD'
            if sim_time - user.last_request >= sim.minimum_period_to_request
                l_best = calculateBestQuality(video(user.video_id), user);
                
                if l_best == user.requested_quality
                    % maintain quality
                elseif l_best > user.requested_quality
                    if user.buffer > sim.qaad.marginal_buffer_length * 1000
                        user.requested_quality = user.requested_quality + 1;
                    end
                else
                    k = 0;
                    
                    while true
                        t = (user.buffer - sim.qaad.minimal_buffer_length*1000) / (1 - (user.estimated_bandwidth / (video(user.video_id).representations_bitrate(user.requested_quality - k)*10^6)));
                        
                        %%%%
                        if t < 0 && (user.buffer > sim.qaad.minimal_buffer_length * 1000)
                            k = k + 1;
                            break
                        end

                        n = (t * user.estimated_bandwidth) / (video(user.video_id).segment_length * video(user.video_id).representations_bitrate(user.requested_quality - k)*10^6);
                        k = k + 1;
                        
                        %if n < 1 && k < user.requested_quality - 1 
                        if n < 1 && k < user.requested_quality
                            continue
                        else
                            break
                        end
                    end
                    
                    %user.requested_quality = user.requested_quality - k;
                    user.requested_quality = user.requested_quality - (k-1);
                end
                
                user.last_request = sim_time;
                user.requesting_data = true;
                user.served_qualities = [user.served_qualities user.requested_quality];
                
                % calculate requested angle
                [~, best_angle_id] = min(abs(video(user.video_id).default_orientations - user.angle));
                user.requested_direction = [user.requested_direction video(user.video_id).default_orientations(best_angle_id)];
                
                user.nr_target_bits = video(user.video_id).nr_bits_per_segment(user.requested_quality);
            end
        case 'QDASH'
            if sim_time - user.last_request >= sim.minimum_period_to_request
                l_best = calculateBestQuality(video(user.video_id), user);
                
                if l_best >= user.requested_quality
                    user.requested_quality = l_best;
                elseif l_best < user.requested_quality
                    if l_best < user.requested_quality - 1
                        t = user.buffer / ((1-user.estimated_bandwidth)/video(user.video_id).representations_bitrate(l_best + 1)*10^6);
                        n = (t * user.estimated_bandwidth)/(video(user.video_id).segment_length * video(user.video_id).representations_bitrate(l_best + 1)*10^6);
                        
                        if n >= 1
                            user.requested_quality = l_best + 1;
                        else
                           user.requested_quality = l_best;
                        end
                    end
                else
                    user.requested_quality = l_best;
                end
                
                user.last_request = sim_time;
                user.requesting_data = true;
                user.served_qualities = [user.served_qualities user.requested_quality];
                
                % calculate requested angle
                [~, best_angle_id] = min(abs(video(user.video_id).default_orientations - user.angle));
                user.requested_direction = [user.requested_direction video(user.video_id).default_orientations(best_angle_id)];
                
                user.nr_target_bits = video(user.video_id).nr_bits_per_segment(user.requested_quality);
            end 
        otherwise
            error('Rate adaptation algorithm is not valid.');
    end
end

% Request better quality if not already at max
function user = requestBetterVideoQuality(user)
    if user.requested_quality < 15
        user.requested_quality = user.requested_quality + 1;
    end
end

% Request worse video quality if not already at min
function user = requestWorseVideoQuality(user)
    if user.requested_quality > 1
        user.requested_quality = user.requested_quality - 1;
    end
end

% Calculate the best possible quality based on estimated bandwidth
function l_best = calculateBestQuality(video, user)
    l_best = 1;

    for representation_id = 1:length(video.representations_bitrate)
        if user.estimated_bandwidth >= video.representations_bitrate(representation_id) * 10^6
            l_best = representation_id;
        else
            break
        end
    end
end

% Update CQI if there is a new CQI event
function user = updateCQIs(sim_time, user)
    for user_id = 1:length(user)
        if user(user_id).cqi_events.event_id > length(user(user_id).cqi_events.time)
            return
        end

        % check for a CQI update only every 5 ms
        if mod(sim_time-1, 5) == 0
             if user(user_id).cqi_events.time(user(user_id).cqi_events.event_id) == sim_time - 1
                user(user_id).cqi = user(user_id).cqi_events.cqi(user(user_id).cqi_events.event_id);
                user(user_id).cqi_events.event_id = user(user_id).cqi_events.event_id + 1;
             end
        end
    end
end

%Allocate RBs
function user = allocateRBs(sim_time, sim, user)
    switch sim.scheduling_mode
        case 'PF1'
            achievable_throughput_1tti = zeros(1,length(user));
            metric = zeros(1,length(user));
            
            if sim_time == 1
                for rb_id = 1:sim.nr_rb_per_tti
                    [~,user_w_less_rbs] = min(extractfield(user, 'nr_allocated_rbs')); % get user with less RBs
                    user(user_w_less_rbs).nr_allocated_rbs = user(user_w_less_rbs).nr_allocated_rbs + 1; % allocate 1 RB to that user
                end
            else
                for user_id = 1:length(user)
                    user(user_id).nr_allocated_rbs = 0;
                    achievable_throughput_1tti(user_id) = sim.nr_rb_per_tti * 1000 * sim.nr_bits_per_rb_per_cqi(user(user_id).cqi);
                    metric(user_id) = achievable_throughput_1tti(user_id)/user(user_id).avg_throughput;
                end
                
                [~,user_ids_w_highest_metric_decrescent] = maxk(metric, length(user));
                
                % allocate RBs to user with highest metric that is requesting data
                for id = 1:length(user)
                    if user(user_ids_w_highest_metric_decrescent(id)).requesting_data == true
                        user(user_ids_w_highest_metric_decrescent(id)).nr_allocated_rbs = sim.nr_rb_per_tti;
                        break
                    end
                end
            end
        case 'PF2'
            % instead of making the calculations for all users and only check if they are valid after, first get valid users, then make calculations
            % reset allocated rbs
            for user_id = 1:length(user)
                user(user_id).nr_allocated_rbs = 0;
            end
            
            % get ids of users that verify the requesting data and latency requirements
            users_requesting_data = [user.requesting_data];
            users_latency = sim_time - [user.last_request] >= sim.latency;
            valid_users = and(users_requesting_data, users_latency);
            valid_users_ids = find(valid_users == true);
            
            achievable_throughput_1rb = zeros(1,length(valid_users_ids));
            metric = zeros(1,length(valid_users_ids));
            
            for id = 1:length(valid_users_ids)
                achievable_throughput_1rb(id) = 1000 * sim.nr_bits_per_rb_per_cqi(user(valid_users_ids(id)).cqi);
                metric(id) = achievable_throughput_1rb(id)/(user(valid_users_ids(id)).avg_throughput);
            end
            
            if ~isempty(valid_users_ids)
                for rb_id = 1:sim.nr_rb_per_tti
                    [~,valid_user_id_w_highest_metric] = maxk(metric, 1);
                    
                    % break cycle if there are no users left needing rbs
                    if isempty(valid_user_id_w_highest_metric)
                        break
                    end
                    
                    user_id = valid_users_ids(valid_user_id_w_highest_metric);
                    
                    user(user_id).nr_allocated_rbs = user(user_id).nr_allocated_rbs + 1;
                    metric(valid_user_id_w_highest_metric) = achievable_throughput_1rb(valid_user_id_w_highest_metric)/(user(user_id).avg_throughput + (user(user_id).nr_allocated_rbs * achievable_throughput_1rb(valid_user_id_w_highest_metric)));
                    
                    % if the received bits + the bits the user will receive with these rbs >= requested bits, stop allocating rbs to this user
                    if user(user_id).nr_received_bits + user(user_id).nr_allocated_rbs * sim.nr_bits_per_rb_per_cqi(user(user_id).cqi) >= user(user_id).nr_target_bits
                        valid_users_ids(valid_user_id_w_highest_metric) = [];
                        metric(valid_user_id_w_highest_metric) = [];
                    end
                end
            end
        case 'RR'
            for user_id = 1:length(user)
                user(user_id).nr_allocated_rbs = 0;
            end
            
            for rb_id = 1:sim.nr_rb_per_tti
                [~,users_w_less_rbs] = mink(extractfield(user, 'nr_allocated_rbs'), length(user)); % get user with less RBs
                
                for id = 1:length(user)
                    user_id = users_w_less_rbs(id);
                    
                    if user(user_id).requesting_data == true && sim_time - user(user_id).last_request >= sim.latency
                        user(user_id).nr_allocated_rbs = user(user_id).nr_allocated_rbs + 1; % allocate 1 RB to that user
                        break
                    end
                end
            end
        case'BET'
            achievable_throughput_1rb = zeros(1,length(user));
            metric = zeros(1,length(user));
    
            for user_id = 1:length(user)
                user(user_id).nr_allocated_rbs = 0;
                metric(user_id) = 1/user(user_id).avg_throughput;
            end

            for rb_id = 1:sim.nr_rb_per_tti
                [~,user_ids_w_highest_metric_decrescent] = maxk(metric, length(user));

                % allocate RBs to user with highest metric that is requesting data
                for id = 1:length(user)
                    user_id = user_ids_w_highest_metric_decrescent(id);

                    if user(user_id).requesting_data == true && sim_time - user(user_id).last_request >= sim.latency
                        user(user_id).nr_allocated_rbs = user(user_id).nr_allocated_rbs + 1;
                        metric(user_id) = 1/((user(user_id).avg_throughput*(sim_time-1) + (user(user_id).nr_allocated_rbs * achievable_throughput_1rb(user_id)))/sim_time);
                        break
                    end
                end
            end
        otherwise
            error('Scheduling algorithm is not valid.');
    end
end

% Write final QoE
function final_qoe = writeFinalQoE(user, final_qoe)
    for user_id = 1:length(user)
        final_qoe(user_id) = user(user_id).qoe;
    end
end

% Write final QoE adjusted by the trajectory
function final_qoe_adjusted = writeFinalQoEAdjusted(user, final_qoe_adjusted)
    for user_id = 1:length(user)
        final_qoe_adjusted(user_id) = user(user_id).adjusted_qoe;
    end
end

% Plot graph with results
function plotGraphs(aux_results)
    figure;
    subplot(4,2,1);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).qoe);
        plot(aux_results(i).adjusted_qoe,'--');
    end
    title('QoE');
    
    subplot(4,2,2);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).buffer);
    end
    title('Buffer');
    
    subplot(4,2,3);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).throughput);
    end
    title('Throughput');
    
    subplot(4,2,4);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).cqi);
    end
    title('CQI');
    
    subplot(4,2,5);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).estimated_bandwidth);
    end
    title('Estimated bandwidth');
    
    subplot(4,2,6);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).requested_quality);
    end
    title('Requested quality');
    
    subplot(4,2,7);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).avg_throughput);
    end
    title('Average throughput');
    
    subplot(4,2,8);
    hold on
    for i = 1:length(aux_results)
        plot(aux_results(i).requested_bitrate);
    end
    title('Requested bitrate');
end