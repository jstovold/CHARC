
%% struct('individual', individual, 'config', config, 'wins', 0, 'visits', 0, 'childrenRemaining', 1:individual.nodes);
    function [selected_idx] = selection(t, node_idx, c)
        
%         this_node           = t.getvalue(node_idx);
        remainingChildren   = t.getRemainingChildren(node_idx);

        % are there any unexpanded child nodes?
        if length(remainingChildren) > 0         % || length(t.getchildren(node_idx)) == 0
            selected_idx  = node_idx;

        else
            % if not, UCB1:
            
            child_idxs  = t.getchildren(node_idx);
            N           = t.getVisits(node_idx); 
            best_idx    = 0;
            best_ucb    = -inf;
            
            for s = child_idxs

                this_w    = t.getWins(s);
                this_n    = t.getVisits(s);
                
                this_ucb  = this_w / this_n + c * sqrt( log(N) / this_n);
                
                if this_ucb > best_ucb
                    best_ucb = this_ucb;
                    best_idx = s;
                end
            end
            if best_idx == 0
                best_idx = randsample(child_idxs, 1);
            end
            [selected_idx] = selection(t, best_idx, c); 
        end
        
    end
    
   