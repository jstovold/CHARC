%%
    function [t] = backprop(t, node_idx, win)
        
%         this_node        = t.getvalue(node_idx);
        wins   = t.getWins(node_idx) + win;
        visits = t.getVisits(node_idx) + 1;
        t.setWins(node_idx, wins);
        t.setVisits(node_idx, visits);
        
        if node_idx ~= 1
            
            parent_idx = t.getparent(node_idx);
            t = backprop(t, parent_idx, win);
            
        end
        
    end


    