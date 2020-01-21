 %% expansion
    % struct('individual', individual, 'config', config, 'wins', 0, 'visits', 0, 'childrenRemaining', 1:individual.nodes);
    % branching factor = num.neurons 
    % add new child nodes to the tree
    
    function [t, new_idx] = expansion(t, node_idx, base_ind)
%         this_node = t.getvalue(node_idx);


        W           = t.getWeights(node_idx);
        Win         = t.getInputs(node_idx);
        Wout        = t.getOutputs(node_idx);
        knockouts   = t.getKnockouts(node_idx);
        children    = t.getRemainingChildren(node_idx);

        if length(children) > 0
            if length(children) == 1
                neuron = children(1);
            else
                neuron = randsample(children, 1);
            end
            
            [W, Win, Wout, knockouts] = knockout_neuron(W, Win, Wout, neuron, knockouts);
%             [ind_prime, conf_prime] = knockout_neuron(ind, neuron, config);
            children(children == neuron) = []; % remove child from remaining children
            
            t.setRemainingChildren(node_idx, children);

            newNodeVal = struct('W', W, ...
                                'input_weights', Win, ...
                                'output_weights', Wout, ...
                                'knockouts', knockouts, ...
                                'wins', 0, ...
                                'visits', 0, ...
                                'childrenRemaining', (setdiff(1:base_ind.nodes, knockouts)));
                            
            new_idx = t.addnode(node_idx, newNodeVal);
            
        end
    end
    
    
    
    