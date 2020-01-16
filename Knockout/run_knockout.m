%% 

function [diff, metrics_a, metrics_b] = run_knockout(individual, config)

    neuron_num = findSynchronisedNeurons(individual, config);
    [ind_prime, config_prime] = knockout_neuron(individual, neuron_num, config);
    [diff, metrics_a, metrics_b] = compareRes(individual, ind_prime, config, config_prime);
    
