%%
function neuron_num = findSynchronisedNeurons(individual, config)

    % run the reservoir with the test sequence and record the output
    collectedStates = collectRoRStates(individual, config.test_input_sequence, config);
    
    % see which neurons are behaving similarly...
    corrs = corrcoef(collectedStates);
    
    % get rid of autocorrelations...
    corrs(corrs == 1) = 0;
    
    % what's the most highly non-auto-correlated correlation
    max_correlation = max(corrs,[],'all');
    
    % and which neuron is this? (there are always two, because it's a correlation between two
    % neurons, but we just pick the first which is returned.
    [row, ~] = find(corrs == max_correlation);
    
    neuron_num = min(row);
    if neuron_num > config.num_nodes
        neuron_num = 0;
    end





