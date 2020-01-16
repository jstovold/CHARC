%% 
function [difference_vec, metrics_a, metrics_b] = compareRes(individual, individual_prime, config, config_prime)

    metrics_a = getMetrics(individual, config);
    metrics_b = getMetrics(individual_prime, config_prime);
    

    % this needs to be extended to also compare their performance on the original task they were
    % evolved for...
    
    test_a = config.testFcn(individual, config);
    test_b = config.testFcn(individual_prime, config_prime);
    
    metrics_a = [metrics_a test_a.test_error];
    metrics_b = [metrics_b test_b.test_error];
    
    difference_vec = metrics_a - metrics_b;

    