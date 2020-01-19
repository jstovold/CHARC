    %% rollout
    % randomly remove neurons until num.neurons = pruneTargetSize
    % evaluate using metricsFcn(.) {metrics(n=10) <= metrics(n = num.neurons)
    
    function [t, win] = rollout(t, node_idx, pruneTargetSize, metricsFcn, baseline, config)
        
%         starting_node = t.getvalue(node_idx);
        starting_ind  = t.getIndividual(node_idx);
        num_neurons   = starting_ind.nodes;
        
        ind  = starting_ind;
        conf = config;
%         conf = t.getConfig(node_idx);
        while num_neurons > pruneTargetSize
            neuron = randsample(1:ind.nodes, 1);
            [ind, conf] = knockout_neuron(ind, neuron, conf);
            num_neurons = num_neurons - 1;
        end
        
        metrics = metricsFcn(ind, conf);
        if won(metrics, baseline)
            win = 1;
        else
            win = 0; 
        end
    end

    function [win] = won(metrics, baseline)
        win = 0;
        if metrics(1) >= baseline(1)
            if metrics(2) <= baseline(2)
                if metrics(3) >= baseline(3)
                    win = 1;
                end
            end
        end
           
    end
    
