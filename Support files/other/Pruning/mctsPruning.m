%%
function [individual, config, t, node, idx] = mctsPruning(individual, metrics_fcn, prune_target_size, comp_budget, config)


    %% setup
    % node format: 1:current individual, 2:config, 3:wins, 4:visits, 5:{remainingChildren array}
    t = tree({individual, config, 0, 0, {1:individual.nodes}}); 
    
    i = 0;
    while i < comp_budget
        
        [node, idx]     = selection(t, 1, 0.1);
        [t, new_idx]    = expansion(t, idx);
        [t, win]        = rollout(t, new_idx, prune_target_size, metrics_fcn, metrics_fcn(individual, config));
        [t]             = backprop(t, new_idx, win);

        i = i + 1;
        disp(i)
    end
    
    % traverse tree to find the leaf node with the best win/visit ratio (this stage can actually be
    % updated to work based on the metrics at the leaf node, but for now I'm sticking with a pretty
    % standard MCTS approach).
    [node, idx] = action_selection(t, 1);
    node        = node{1};
    individual  = node{1};
    config      = node{2};
    
    %% selection
    % UCB1 = w_i / n_i + c * \sqrt{log(N_i) / n_i}
    
    function [selected_node, selected_idx] = action_selection(t, node_idx)
        this_node           = t.Node(node_idx);
        this_node           = this_node{1};
        % is this a leaf node?
        if t.isleaf(node_idx)
            selected_idx  = node_idx;
            selected_node = t.Node(selected_idx);
        else
            
            child_idxs  = t.getchildren(node_idx);
            N           = this_node{4}; 
            best_idx    = 0;
            best_opt    = -inf;
            
            for s = child_idxs
                this_node = t.Node(s);
                this_node = this_node{1};
                this_w    = this_node{3};
                this_n    = this_node{4};
                
                if (this_w / this_n) > best_opt
                    best_opt = this_w / this_n;
                    best_idx = s;
                end
            end
            if best_idx == 0
                % no children have been visited at this stage
                % pick at random?
                [selected_node, selected_idx] = action_selection(t, randsample(child_idxs, 1));
            else
                [selected_node, selected_idx] = action_selection(t, best_idx); 
            end
        end
    end
    

    function [selected_node, selected_idx] = selection(t, node_idx, c)
        
        this_node           = t.Node(node_idx);
        this_node           = this_node{1};
        remainingChildren   = this_node{5};
%         disp(remainingChildren)
        % are there any unexpanded child nodes?
        if length(remainingChildren{1}) > 0         % || length(t.getchildren(node_idx)) == 0
            selected_idx  = node_idx;
            selected_node = t.Node(selected_idx);
            selected_node = selected_node{1};
%            disp(length(remainingChildren{1}))
        else
            % if not, UCB1:
            
            child_idxs  = t.getchildren(node_idx);
            N           = this_node{4}; 
            best_idx    = 0;
            best_ucb    = -inf;
            
            for s = child_idxs
                this_node = t.Node(s);
                this_node = this_node{1};
                this_w    = this_node{3};
                this_n    = this_node{4};
                
                this_ucb  = this_w / this_n + c * sqrt( log(N) / this_n);
                
                if this_ucb > best_ucb
                    best_ucb = this_ucb;
                    best_idx = s;
                end
            end
            if best_idx == 0
                best_idx = randsample(child_idxs, 1);
            end
            [selected_node, selected_idx] = selection(t, best_idx, c); 
        end
        
    end
    
    %% expansion
    % branching factor = num.neurons 
    % add new child nodes to the tree
    
    function [t, new_idx] = expansion(t, node_idx)
        this_node = t.Node(node_idx);
        this_node = this_node{1};
        ind       = this_node{1};
        conf      = this_node{2};
        children  = this_node{5};
%         disp(children)
        if length(children{1}) > 0
            if length(children{1}) == 1
                neuron = children{1}(1);
            else
                neuron = randsample(children{1}, 1);
            end
            [ind_prime, conf_prime] = knockout_neuron(ind, neuron, conf);
            children{1}(children{1} == neuron) = []; % remove child from remaining children
            this_node{5} = children;
            t = t.set(node_idx, this_node);
            [t, new_idx] = t.addnode(node_idx, {ind_prime, conf_prime, 0, 0, {1:ind_prime.nodes}});
        end
    end
    
    
    
    
    %% rollout
    % randomly remove neurons until num.neurons = pruneTargetSize
    % evaluate using metricsFcn(.) {metrics(n=10) <= metrics(n = num.neurons)
    
    function [t, win] = rollout(t, node_idx, pruneTargetSize, metricsFcn, baseline)
        
        starting_node = t.Node(node_idx);
        starting_node = starting_node{1};
        starting_ind  = starting_node{1};
        num_neurons   = starting_ind.nodes;
        
        ind  = starting_ind;
        conf = starting_node{2};
        while num_neurons > pruneTargetSize
            neuron = randsample(1:ind.nodes, 1);
            [ind, conf] = knockout_neuron(ind, neuron, conf);
            num_neurons = num_neurons - 1;
        end
        
        metrics = metricsFcn(ind, conf);
        if wins(metrics, baseline)
            win = 1;
        else
            win = 0; 
        end
    end

    function [win] = wins(metrics, baseline)
        win = 0;
        if metrics(1) >= baseline(1)
            if metrics(2) <= baseline(2)
                if metrics(3) >= baseline(3)
                    win = 1;
                end
            end
        end
           
    end
    
    
    %% backpropagation
    % update tree


    function [t] = backprop(t, node_idx, win)
        
        this_node    = t.Node(node_idx);
        this_node    = this_node{1};
        this_node{3} = this_node{3} + win;
        this_node{4} = this_node{4} + 1;
        t = t.set(node_idx, this_node);
        
        if node_idx ~= 1
            
            parent_idx = t.getparent(node_idx);

            t = backprop(t, parent_idx, win);
            
        end
        
    end





end




