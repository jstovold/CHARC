%%
function [selected_node, selected_idx] = action_selection(t, node_idx)
%     this_node           = t.getvalue(node_idx);
    
    % is this a leaf node?
    if t.isleaf(node_idx)
        selected_idx  = node_idx;
        selected_node = t.getnode(selected_idx);
    else

        child_idxs  = t.getchildren(node_idx);
        best_idx    = 0;
        best_w  = -inf;

        for s = child_idxs
%             this_node = t.getvalue(s);
            this_w    = t.getWins(s);
            this_n    = t.getVisits(s);

            if this_w > best_w

                best_w = this_w;
                best_idx = s;
            end

        end
        if best_idx == 0
            % no children have been visited at this stage
            % pick at random?
            disp('child');
            disp(length(child_idxs))
            [selected_node, selected_idx] = action_selection(t, randsample(child_idxs, 1));
        else
            disp(best_idx)
            disp(best_w)
            [selected_node, selected_idx] = action_selection(t, best_idx); 
        end

    end
end
    