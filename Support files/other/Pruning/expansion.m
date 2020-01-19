 %% expansion
    % struct('individual', individual, 'config', config, 'wins', 0, 'visits', 0, 'childrenRemaining', 1:individual.nodes);
    % branching factor = num.neurons 
    % add new child nodes to the tree
    
    function [t, new_idx] = expansion(t, node_idx, config)
%         this_node = t.getvalue(node_idx);

        ind       = t.getIndividual(node_idx);
%         conf      = t.getConfig(node_idx);
        children  = t.getRemainingChildren(node_idx);

        if length(children) > 0
            if length(children) == 1
                neuron = children(1);
            else
                neuron = randsample(children, 1);
            end
            [ind_prime, conf_prime] = knockout_neuron(ind, neuron, config);
            children(children == neuron) = []; % remove child from remaining children
            
            t.setRemainingChildren(node_idx, children);
            
%             this_node.childrenRemaining = children;
%             t.updateValue(node_idx, this_node);
            newNodeVal = struct('individual', ind_prime, 'knockouts', conf_prime.knocked_out_neurons, 'wins', 0, 'visits', 0, 'childrenRemaining', 1:ind_prime.nodes);
            new_idx = t.addnode(node_idx, newNodeVal);
            
        end
    end
    
    
    
    