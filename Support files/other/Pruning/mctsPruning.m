%%
function [individual, config, t, node, idx, total_wins] = mctsPruning(individual, metrics_fcn, prune_target_size, comp_budget, config)

    % node format: 1:current individual, 2:config, 3:wins, 4:visits, 5:{remainingChildren array}
    
    if ~isfield(config, 'knocked_out_neurons')
        config.knocked_out_neurons = [];
    end
    rootVal = struct('individual', individual, 'knockouts', [], 'wins', 0, 'visits', 0, 'childrenRemaining', 1:individual.nodes);
    t = jhs_tree(comp_budget);
    rootID = t.addroot(rootVal);
    
%     t = tree({individual, config, 0, 0, {1:individual.nodes}}); 
    baseline = metrics_fcn(individual, config);
    i = 0;
    while i < comp_budget
        
        [idx]           = selection(t, rootID, 0.1);
        [t, new_idx]    = expansion(t, idx, config);
        [t, win]        = rollout(t, new_idx, prune_target_size, metrics_fcn, baseline, config);
        [t]             = backprop(t, new_idx, win);

        i = i + 1;
        disp(i)
    end
    
    % traverse tree to find the leaf node with the best win/visit ratio (this stage can actually be
    % updated to work based on the metrics at the leaf node, but for now I'm sticking with a pretty
    % standard MCTS approach).
    [node, idx]                     = action_selection(t, 1);
    node                            = node.value;
    individual                      = node.individual;
    config.knocked_out_neurons      = node.knockouts;
    root                            = t.getvalue(1);
    total_wins                      = root.wins;
 



end




