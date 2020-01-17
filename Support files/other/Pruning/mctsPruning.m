%%
function [individual] = mctsPruning(individual, metrics_fcn, prune_target_size, comp_budget, config)

    c = 0.1;
    % this is just to shut up the IDE:
    a = metrics_fcn(prune_target_size, comp_budget, config);
    disp(a)

    %% setup
    t = tree({individual, config, 0, 0, 0}); % root node: current individual, config, wins, visits, UCB value
    layer1_nodes = zeros(individual.nodes, 1);
    for n = 1:individual.nodes
        [individual_prime, config_prime] = knockout_neuron(individual, n, config);
        [t, layer1_nodes(n)] = t.addnode(1, {individual_prime, config_prime, 0, 0, 0}); 
    end
    
    %% selection
    % UCB1 = w_i / n_i + c * \sqrt{log(N_i) / n_i}
    
    
    function [selected_node, selected_idx] = selection(t, node_idx, c)
        % are there any unexpanded child nodes?
        this_node  = t.Node(node_idx);
        ind        = this_node{1};
        
        if length(t.getchildren(node_idx)) < ind.nodes
            selected_idx  = node_idx;
            selected_node = t.Node(selected_idx);
        else
            % if not, UCB1:
            best_idx    = 0;
            best_ucb    = -inf;
            child_nodes = t.getchildren(node_idx);
            
            % JHS continue here
            % this still needs stripping apart and re-writing so that the traversal functions
            % correctly.
            
            for s = 1:child_nodes
                parent_node = t.get(t.getparent(s));
                N = parent_node{4};
                this_node = t.Node(s);
                w_i = this_node{3};
                n_i = this_node{4};

                % update ucb value in the tree
                new_ucb = w_i / n_i + c * sqrt( log(N) / n_i);
                this_node{5} = new_ucb;
                t = t.set(s, this_node);

                if new_ucb > best_ucb
                    best_ucb = new_ucb;
                    best_idx = s;
                end
            end

            % best_idx should contain the index of our selected node
            selected_idx  = best_idx;
            selected_node = t.Node(selected_idx);
        end
    end
    
    %% expansion
    % branching factor = num.neurons 
    % add new child nodes to the tree
    layer2_nodes = zeros(selected_ind.nodes, 1);
    for n = 1:selected_ind.nodes
        [individual_prime, config_prime] = knockout_neuron(selected_ind, n, config);
        [t, layer2_nodes(n)] = t.addnode(selected_idx, {individual_prime, config_prime, 0, 0, 0}); 
    end
    
    
    
    %% rollout
    % randomly remove neurons until num.neurons = pruneTargetSize
    % evaluate using metricsFcn(.) {metrics(n=10) <= metrics(n = num.neurons)
    
    %% backpropagation
    % update tree













