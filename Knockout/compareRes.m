%% 
function difference_vec = compareRes(individual, individual_prime, config)

    metrics_a = getMetrics(individual, config);
    metrics_b = getMetrics(individual_prime, config);
    
    difference_vec = metrics_a - metrics_b;

    
    
    % this needs to be extended to also compare their performance on the original task they were
    % evolved for...
    
    
    
    
    