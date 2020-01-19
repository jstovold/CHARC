
%%
classdef jhs_tree < handle
   properties 
      treeStructure(1,:) {mustBeNumeric}    % vector of parentIDs, stores the structure of the tree
%       valueStructure(1,:)                   % vector of nodes, stores the underlying values for the tree
      idStructure(1,:)
      individualStructure(1,:)
      knockoutsStructure(1, :)
      winsStructure(1,:)
      visitsStructure(1,:)
      remainingChildrenStructure(1,:)
      childrenIDsStructure(1,:)
      nextEmptyIdx
   end
   methods
       function obj = jhs_tree(size)
           obj.nextEmptyIdx = 1;
           obj.treeStructure  = zeros(1, size + 1);
           
           
           %obj.valueStructure = %repmat(struct('parent', 0, 'ID', 0, 'value', {0}, 'childrenIDs', []), size, 1 );
       end
       function rootID = addroot(obj, val)
            rootID = obj.nextEmptyIdx;
            
            obj.treeStructure(rootID) = 0;
%             obj.valueStructure = repmat(struct('parent', 0, 'ID', rootID, 'value', val, 'childrenIDs', []), length(obj.treeStructure), 1 );
            obj.idStructure                 = repmat([0], length(obj.treeStructure), 1);
            obj.individualStructure         = repmat([val.individual], length(obj.treeStructure), 1);
            obj.knockoutsStructure          = repmat({[]}, length(obj.treeStructure), 1);
            obj.winsStructure               = repmat([val.wins], length(obj.treeStructure), 1);
            obj.visitsStructure             = repmat([val.visits], length(obj.treeStructure), 1);
            obj.remainingChildrenStructure  = repmat({val.childrenRemaining}, length(obj.treeStructure), 1);
            obj.childrenIDsStructure        = repmat({[]}, length(obj.treeStructure), 1);
            
            obj.idStructure(rootID)         = rootID;
%             obj.valueStructure(rootID) = struct('parent', 0, 'ID', rootID, 'value', val, 'childrenIDs', []);
            obj.nextEmptyIdx = obj.nextEmptyIdx + 1;
           
       end
       function nodeID = addnode(obj, parentID, val)
           if parentID == 0
               nodeID = obj.addroot(val);
           else
               nodeID = obj.nextEmptyIdx;
               if nodeID > length(obj.treeStructure)
                   error('Not enough space to add new node');
               end
               obj.treeStructure(nodeID) = parentID;
%                obj.valueStructure(nodeID) = struct('parent', parentID, 'ID', nodeID, 'value', val, 'childrenIDs', []);
               
               obj.idStructure(nodeID)                  = nodeID;
               obj.individualStructure(nodeID)          = val.individual;
               obj.knockoutsStructure(nodeID)           = {val.knockouts};
               obj.winsStructure(nodeID)                = val.wins;
               obj.visitsStructure(nodeID)              = val.visits;
               obj.remainingChildrenStructure(nodeID)   = {val.childrenRemaining};
%                obj.childrenIDsStructure(nodeID)         = [];
               
               childIDs                                 = obj.childrenIDsStructure(parentID);
               obj.childrenIDsStructure(parentID)       = {[childIDs{1}, nodeID]};
               obj.nextEmptyIdx = obj.nextEmptyIdx + 1;
           end
       end
       function val = getvalue(obj,idx)
           val = obj.getnode(idx).value;
       end
       function node = getnode(obj,idx)
%            node = obj.valueStructure(idx);
           valueStruct = struct('individual', obj.individualStructure(idx), ...
                                'knockouts', obj.knockoutsStructure(idx), ...
                                'wins', obj.winsStructure(idx), ...
                                'visits', obj.visitsStructure(idx), ...
                                'childrenRemaining', obj.remainingChildrenStructure(idx));
           node = struct('parent', obj.treeStructure(idx), 'ID', idx, 'value', valueStruct, 'childrenIDs', obj.childrenIDsStructure(idx));
       end
       
       function individual = getIndividual(obj, idx)
           individual = obj.individualStructure(idx);
       end
       
       function knockouts = getKnockouts(obj, idx)
           knockouts = obj.knockoutsStructure(idx);
           knockouts = knockouts{1};

       end
       
       function wins = getWins(obj, idx)
           wins = obj.winsStructure(idx);
       end
       function visits = getVisits(obj, idx)
           visits = obj.visitsStructure(idx);
       end
       function remainingChildren = getRemainingChildren(obj, idx)
           remainingChildren = obj.remainingChildrenStructure(idx);
           remainingChildren = remainingChildren{1};
       end
       
         
       function setIndividual(obj, idx, val)
           obj.individualStructure(idx) = val;
       end
       
       function setKnockouts(obj, idx, val)
           obj.knockoutsStructure(idx) = {val};
       end
       
       function setWins(obj, idx, val)
           obj.winsStructure(idx) = val;
       end
       function setVisits(obj, idx, val)
           obj.visitsStructure(idx) = val;
       end
       function setRemainingChildren(obj, idx, val)
           obj.remainingChildrenStructure(idx) = {val};
       end
       
       
       
       function parent = getparent(obj, idx)
           parent = obj.treeStructure(idx);
       end
       function [children] = getchildren(obj, idx)
           children = obj.childrenIDsStructure(idx);
           children = children{1};
%            children = obj.valueStructure(idx).childrenIDs;
       end
%        function updateValue(obj, idx, newval)
%            nodeVal = obj.valueStructure(idx);
%            nodeVal.value = newval;
%            obj.valueStructure(idx) = nodeVal;
%        end

       function leaf = isleaf(obj, idx)
           leaf = isempty(obj.getchildren(idx));
       end
       
   end
end

% node format: {struct: parentID | ID | value | [childrenIDs] }
% tree format: {vector: parentID | parentID | parentID ... }


