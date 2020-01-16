%% 

function output = run_knockout(individual, config)

    neuron_num = findSynchronisedNeurons(individual, config);
    [ind_prime, config_prime] = knockout_neuron(individual, neuron_num, config);
    diff = compareRes(individual, ind_prime, config);
    output = diff;
