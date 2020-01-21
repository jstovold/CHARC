    %% rollout
    % randomly remove neurons until num.neurons = pruneTargetSize
    % evaluate using metricsFcn(.) {metrics(n=10) <= metrics(n = num.neurons)
    
    function [t, win] = rollout(t, node_idx, pruneTargetSize, metricsFcn, baseline, config)
        
%         starting_node = t.getvalue(node_idx);
        [starting_ind, knockouts] = t.getIndividual(node_idx);
%         num_neurons = starting_ind.nodes;
        
        ind  = starting_ind;
        conf = config;
        W    = ind.W;
        Win  = ind.input_weights;
        Wout = ind.output_weights{1};
        
        neurons_in_ind = setdiff(1:ind.nodes, knockouts{1});
        num_neurons = length(neurons_in_ind);
        while num_neurons > pruneTargetSize
            neurons_in_ind = setdiff(1:ind.nodes, knockouts{1});
            neuron = randsample(neurons_in_ind, 1);
            [W, Win, Wout, knockouts] = knockout_neuron(W, Win, Wout, neuron, knockouts);
            num_neurons = num_neurons - 1;
        end
        ind_p = t.buildIndividual(ind, W, Win, Wout, knockouts);
%         disp(ind_p)
        metrics = metricsFcn(ind_p, conf);
        if won(metrics, baseline)
            win = 1;
        else
            win = 0; 
        end
    end

    function [win] = won_old(metrics, baseline)
        win = 0;
        if metrics(1) >= baseline(1)
            if metrics(2) <= baseline(2)
                if metrics(3) >= baseline(3)
                    win = 1;
                end
            end
        end
           
    end
    
    function [win] = won(metrics, baseline)
        win = 0;
        if metrics(1) == baseline(1)
            if metrics(2) == baseline(2)
                if abs(metrics(3) - baseline(3)) <= 0.5
                    win = 1;
                end
            end
        end
    
    end
