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
    total_start = tic;
    reverseLen = 0;
    msg = '';
    while i < comp_budget
        a = tic;
        [idx]           = selection(t, rootID, 0.1);
        [t, new_idx]    = expansion(t, idx, individual);
        [t, win]        = rollout(t, new_idx, prune_target_size, metrics_fcn, baseline, config);
        [t]             = backprop(t, new_idx, win);

        i = i + 1;
        
        fprintf(strcat(repmat('\b', 1, reverseLen)));
        reverseLen = fprintf("Iter: %i \t Time: %.3fs ", i, toc(a));
        
        
        
    end
    total_time = toc(total_start);
    
    % traverse tree to find the leaf node with the best win/visit ratio (this stage can actually be
    % updated to work based on the metrics at the leaf node, but for now I'm sticking with a pretty
    % standard MCTS approach).
    [node, idx]                     = action_selection(t, 1);
    node                            = node.value;
    individual                      = t.getIndividual(idx);
    config.knocked_out_neurons      = node.knockouts;
    root                            = t.getvalue(1);
    total_wins                      = root.wins;
 
    fprintf('\r\nTotal time for %i iterations: %fmins \n\r', i, total_time/60);


end




