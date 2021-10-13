clear
rng shuffle
rand_stream = rng;

%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%
% parameters(1).latency = 10; % request-response latency - in ms (fixed value)
% parameters(1).scenario = 'monolithic'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(1).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(1).buffer_request_threshold = 6000; % buffer level that triggers a new request - in ms
% parameters(1).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(2).latency = 1; % request-response latency - in ms (fixed value)
% parameters(2).scenario = 'monolithic'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(2).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(2).buffer_request_threshold = 6000; % buffer level that triggers a new request - in ms
% parameters(2).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(3).latency = 10; % request-response latency - in ms (fixed value)
% parameters(3).scenario = 'tiles'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(3).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(3).buffer_request_threshold = [50 100 300 500 700 900 1000 2000 3000 4000 6000]; % buffer level that triggers a new request - in ms
% parameters(3).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(4).latency = 1; % request-response latency - in ms (fixed value)
% parameters(4).scenario = 'tiles'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(4).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(4).buffer_request_threshold = [50 100 300 500 700 900 1000 2000 3000 4000 6000]; % buffer level that triggers a new request - in ms
% parameters(4).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(5).latency = 10; % request-response latency - in ms (fixed value)
% parameters(5).scenario = 'tiles_partial'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(5).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(5).buffer_request_threshold = [50 100 300 500 700 900 1000 6000]; % buffer level that triggers a new request - in ms
% parameters(5).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(6).latency = 1; % request-response latency - in ms (fixed value)
% parameters(6).scenario = 'tiles_partial'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(6).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(6).buffer_request_threshold = [50 100 300 500 700 900 1000 6000]; % buffer level that triggers a new request - in ms
% parameters(6).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(7).latency = 10; % request-response latency - in ms (fixed value)
% parameters(7).scenario = 'viewport_only'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(7).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(7).buffer_request_threshold = [1 2 3 4 5 10 15 20 30 40]; % buffer level that triggers a new request - in ms
% parameters(7).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(8).latency = 1; % request-response latency - in ms (fixed value)
% parameters(8).scenario = 'viewport_only'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(8).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
% parameters(8).buffer_request_threshold = [1 2 3 4 5 10 15 20 30 40]; % buffer level that triggers a new request - in ms
% parameters(8).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(9).latency = 10; % request-response latency - in ms (fixed value)
% parameters(9).scenario = 'viewport_only_modified'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(9).nr_users = [1 10 20 30 40 50 60 70 80 90 100 120 140 160]; % number of users
% parameters(9).buffer_request_threshold = [15 20]; % buffer level that triggers a new request - in ms
% parameters(9).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower
% 
% parameters(10).latency = 1; % request-response latency - in ms (fixed value)
% parameters(10).scenario = 'viewport_only_modified'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
% parameters(10).nr_users = [1 10 20 30 40 50 60 70 80 90 100 120 140 160]; % number of users
% parameters(10).buffer_request_threshold = [3 4 5]; % buffer level that triggers a new request - in ms
% parameters(10).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower

parameters(1).latency = 10; % request-response latency - in ms (fixed value)
parameters(1).scenario = 'viewport_only_modified'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
parameters(1).nr_users = [1 10 20 30 40 50 60 70 80 90 100 120 140 160]; % number of users
parameters(1).buffer_request_threshold = [15 20]; % buffer level that triggers a new request - in ms
parameters(1).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower

parameters(2).latency = 1; % request-response latency - in ms (fixed value)
parameters(2).scenario = 'viewport_only_modified'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
parameters(2).nr_users = [1 10 20 30 40 50 60 70 80 90 100 120 140 160]; % number of users
parameters(2).buffer_request_threshold = [3 4 5]; % buffer level that triggers a new request - in ms
parameters(2).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower

parameters(3).latency = 10; % request-response latency - in ms (fixed value)
parameters(3).scenario = 'tiles'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
parameters(3).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
parameters(3).buffer_request_threshold = [2000 3000 4000]; % buffer level that triggers a new request - in ms
parameters(3).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower

parameters(4).latency = 10; % request-response latency - in ms (fixed value)
parameters(4).scenario = 'tiles'; % simulation scenario - monolithic, tiles, tiles_partial, viewport_only
parameters(4).nr_users = [1 10 20 30 40 50 60 70 80 90 100]; % number of users
parameters(4).buffer_request_threshold = [2000 3000 4000]; % buffer level that triggers a new request - in ms
parameters(4).debug_mode = false; % if true, calculate the cumulative QoE and other metric every TTI and plot the graphs, but simulations run slower

%%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%
mc.nr_monte_carlo_sims = 24;

%%%%%%%%%%%%%%%%% RANDOM ALLOCATIONS %%%%%%%%%%%%%%%%%%
allocations = runAllocationAlgorithms(mc, parameters);

%%%%%%%%%%%%%%%%%%% RUN SIMULATIONS %%%%%%%%%%%%%%%%%%%
for sim_id = 1:length(parameters)
    results = runSimulator(allocations, mc, parameters(sim_id));
    save(getFileName(parameters(sim_id)));
end

%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%

function video_allocation = videoAllocation(nr_users)
    video_allocation = randi(3, 1, nr_users); % allocation of one of the videos to each user
end

function trajectory_allocation = trajectoryAllocation(nr_users)
    trajectory_allocation = randperm(48);
    aux = trajectory_allocation;
    nr_remaining_users = nr_users - 48;

    while nr_remaining_users > 0
        trajectory_allocation = [trajectory_allocation aux];
        nr_remaining_users = nr_remaining_users - 48;
    end

    trajectory_allocation = trajectory_allocation(1:nr_users); % allocation of one of the trajectories to each user
end

function cqi_allocation = CQIAllocation(nr_users)
    cqi_allocation = randperm(200);
    aux = cqi_allocation;
    nr_remaining_users = nr_users - 200;

    while nr_remaining_users > 0
        cqi_allocation = [cqi_allocation aux];
        nr_remaining_users = nr_remaining_users - 200;
    end

    cqi_allocation = cqi_allocation(1:nr_users); % allocation of one of the CQIs to each user
end

function initial_buffering_time = initialBufferingTime(nr_users)
    max_time = 200; % mudar tambem na funcao initialBuffering do simulador
    initial_buffering_time = randi(max_time, 1, nr_users);
end

function allocations = runAllocationAlgorithms(mc, parameters)
    max_users = max(parameters(1).nr_users);
    for i = 1:length(parameters)
        if max(parameters(i).nr_users) > max_users
            max_users = max(parameters(i).nr_users);
        end
    end

    allocations.video_allocation = zeros(mc.nr_monte_carlo_sims, max_users);
    allocations.trajectory_allocation = zeros(mc.nr_monte_carlo_sims, max_users);
    allocations.cqi_allocation = zeros(mc.nr_monte_carlo_sims, max_users);
    allocations.initial_buffering_time = zeros(mc.nr_monte_carlo_sims, max_users);

    for id_sim = 1:mc.nr_monte_carlo_sims
        allocations.video_allocation(id_sim,:) = videoAllocation(max_users); % allocation of one of the videos to each user
        allocations.trajectory_allocation(id_sim,:) = trajectoryAllocation(max_users); % allocation of one of the trajectories to each user
        allocations.cqi_allocation(id_sim,:) = CQIAllocation(max_users); % allocation of one of the CQI files to each user
        allocations.initial_buffering_time(id_sim,:) = initialBufferingTime(max_users); % set the time for each user to start initial buffering
    end
end

function results = runSimulator(allocations, mc, parameters)
    nr_monte_carlo_sims = mc.nr_monte_carlo_sims;
    length_buffer_threshold = length(parameters.buffer_request_threshold);
    length_users = length(parameters.nr_users);
    
    %%%%%%%%%%%%%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%
    results.nr_users = NaN;
    results.latency = NaN;
    results.scenario = NaN;
    results.buffer_request_threshold = NaN;
    results.qoe = NaN;
    results.qoe_adjusted = NaN;
    results.global_qoe = NaN;
    results.global_qoe_adjusted = NaN;
    results.satisfaction_rate_3 = NaN;
    results.satisfaction_rate_4 = NaN;
    results.insatisfaction_rate_2 = NaN;
    results.satisfaction_rate_3_adj = NaN;
    results.satisfaction_rate_4_adj = NaN;
    results.insatisfaction_rate_2_adj = NaN;
    results = repmat(results, mc.nr_monte_carlo_sims, length_users, length_buffer_threshold);
    
    %%%%%%%%%%%%%%% MONTE CARLO SIMULATION %%%%%%%%%%%%%%%%
    parfor id_sim = 1:nr_monte_carlo_sims
        for id_buffer_threshold = 1:length_buffer_threshold
            for id_users = 1:length_users
                % write parameters to results
                results(id_sim,id_users,id_buffer_threshold).nr_users = parameters.nr_users(id_users);
                results(id_sim,id_users,id_buffer_threshold).latency = parameters.latency;
                results(id_sim,id_users,id_buffer_threshold).scenario = parameters.scenario;
                results(id_sim,id_users,id_buffer_threshold).buffer_request_threshold = parameters.buffer_request_threshold(id_buffer_threshold);
                results(id_sim,id_users,id_buffer_threshold).qoe = NaN(1, parameters.nr_users(id_users));
                results(id_sim,id_users,id_buffer_threshold).qoe_adjusted = NaN(1, parameters.nr_users(id_users));
                
                fprintf('Simulation %.0f (of %.0f): Nr. users - %.0f; Latency - %.0f; Min. buffer level - %.0f; Scenario - %s\n', id_sim, mc.nr_monte_carlo_sims, parameters.nr_users(id_users), parameters.latency, parameters.buffer_request_threshold(id_buffer_threshold), parameters.scenario);
                
                nr_users = parameters.nr_users(id_users);
                buffer_request_threshold = parameters.buffer_request_threshold(id_buffer_threshold);
                scenario = parameters.scenario;
                latency = parameters.latency;
                debug_mode = parameters.debug_mode;
                
                video_allocation = allocations.video_allocation(id_sim,:);
                trajectory_allocation = allocations.trajectory_allocation(id_sim,:);
                cqi_allocation = allocations.cqi_allocation(id_sim,:);
                initial_buffering_time = allocations.initial_buffering_time(id_sim,:);
                
                % run simulator
                [qoe, qoe_adjusted] = Simulator(nr_users, buffer_request_threshold, scenario, latency, debug_mode, video_allocation, trajectory_allocation, cqi_allocation, initial_buffering_time);

                % write results to results
                results(id_sim,id_users,id_buffer_threshold).qoe = qoe;
                results(id_sim,id_users,id_buffer_threshold).qoe_adjusted = qoe_adjusted;

                % clean results
                results(id_sim,id_users,id_buffer_threshold).qoe(isnan(results(id_sim,id_users,id_buffer_threshold).qoe)) = 0;
                results(id_sim,id_users,id_buffer_threshold).qoe_adjusted(isnan(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted)) = 0;

                % satisfaction calculations
                results(id_sim,id_users,id_buffer_threshold).global_qoe = mean(results(id_sim,id_users,id_buffer_threshold).qoe);
                results(id_sim,id_users,id_buffer_threshold).global_qoe_adjusted = mean(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted);
                results(id_sim,id_users,id_buffer_threshold).satisfaction_rate_3 = length(results(id_sim,id_users,id_buffer_threshold).qoe(results(id_sim,id_users,id_buffer_threshold).qoe >= 3)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
                results(id_sim,id_users,id_buffer_threshold).satisfaction_rate_4 = length(results(id_sim,id_users,id_buffer_threshold).qoe(results(id_sim,id_users,id_buffer_threshold).qoe >= 4)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
                results(id_sim,id_users,id_buffer_threshold).insatisfaction_rate_2 = length(results(id_sim,id_users,id_buffer_threshold).qoe(results(id_sim,id_users,id_buffer_threshold).qoe <= 2)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
                results(id_sim,id_users,id_buffer_threshold).satisfaction_rate_3_adj = length(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted >= 3)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
                results(id_sim,id_users,id_buffer_threshold).satisfaction_rate_4_adj = length(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted >= 4)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
                results(id_sim,id_users,id_buffer_threshold).insatisfaction_rate_2_adj = length(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted(results(id_sim,id_users,id_buffer_threshold).qoe_adjusted <= 2)) / results(id_sim,id_users,id_buffer_threshold).nr_users;
            end
        end
    end
end

function save_file = getFileName(parameters)
    save_file = strcat(parameters.scenario, '_latency_', int2str(parameters.latency), '_bufferlevel');
    for i = 1:length(parameters.buffer_request_threshold)
        save_file = strcat(save_file, '_', int2str(parameters.buffer_request_threshold(i)));
    end
end