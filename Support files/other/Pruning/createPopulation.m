%% 
% creates a population of ESNs (which are currently not trained to any task). for a population 
% which is trained to a specific task, use create_population.m in the knockout folder.

function [population, config, metrics] = createPopulation(num_pop, res_size)



%% Setup
rng(1,'twister');

% type of network to evolve
config.res_type = 'RoR';                    % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [res_size];              % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. 

config = selectReservoirType(config);       % collect function pointers for the selected reservoir type 


%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.preprocess = 1;             % basic preprocessing, e.g. scaling and mean variance
config.dataset = 'NARMA10';        % Task to evolve for

config.prune_iterations = 50;

% get dataset information
[config] = selectDataset(config);

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc. 
[config] = getAdditionalParameters(config);

%% general params
config.gen_print = 25;                       % after 'gen_print' generations print task performance and show any plots
config.start_time = datestr(now, 'HH:MM:SS');
config.save_gen = inf;                       % save data at generation = save_gen
config.pop_size = num_pop;

% type of metrics to apply; if necessary
config.metrics = {'KR','GR','MC'};          % list metrics to apply in cell array: see getVirtualMetrics.m for types of metrics available
config.record_metrics = 0;                  % save metrics


population = config.createFcn(config);
metrics = zeros(num_pop, length(config.metrics));
for i = 1:num_pop
    metrics(i,:) = getMetrics(population(i), config);
end










