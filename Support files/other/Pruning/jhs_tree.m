
%%
classdef jhs_tree < handle
   properties 
      treeStructure(1,:) {mustBeNumeric}    % vector of parentIDs, stores the structure of the tree
%       valueStructure(1,:)                   % vector of nodes, stores the underlying values for the tree
      idStructure(1,:)
%       individualStructure(1,:)
      weightStructure(1,:)
      inputStructure(1,:)
      outputStructure(1,:)
      knockoutsStructure(1, :)
      winsStructure(1,:)
      visitsStructure(1,:)
      remainingChildrenStructure(1,:)
      childrenIDsStructure(1,:)
      individual_base
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
            obj.weightStructure             = repmat(val.individual.W, length(obj.treeStructure), 1);
            obj.inputStructure              = repmat(val.individual.input_weights, length(obj.treeStructure), 1);
            obj.outputStructure             = repmat({val.individual.output_weights}, length(obj.treeStructure), 1);
            obj.individual_base             = val.individual;
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
               obj.weightStructure(nodeID)              = {val.W};
               obj.inputStructure(nodeID)               = {val.input_weights};
               obj.outputStructure(nodeID)              = {val.output_weights};
               obj.knockoutsStructure(nodeID)           = {val.knockouts};
               obj.winsStructure(nodeID)                = val.wins;
               obj.visitsStructure(nodeID)              = val.visits;
               obj.remainingChildrenStructure(nodeID)   = {val.childrenRemaining};
               
               childIDs                                 = obj.childrenIDsStructure(parentID);
               obj.childrenIDsStructure(parentID)       = {[childIDs{1}, nodeID]};
               obj.nextEmptyIdx = obj.nextEmptyIdx + 1;
           end
       end
       function val = getvalue(obj,idx)
           val = obj.getnode(idx).value;
       end
       
       function ind = buildIndividual(obj, base, W, Win, Wout, knockouts)
           
           ind = base;
           ind.W{1,1} = W{1,1};
           ind.input_weights{1,1} = Win{1,1};
           ind.output_weights = Wout;
%            ind.nodes = ind.nodes - length(knockouts);
%            ind.total_units = ind.nodes;
%            
       end
       
       function [ind, knockouts] = getIndividual(obj, idx)
           knockouts = obj.knockoutsStructure(idx);
           ind = obj.buildIndividual(obj.individual_base, ...
                                     obj.weightStructure(idx), ...
                                     obj.inputStructure(idx), ...
                                     obj.outputStructure(idx), ...
                                     knockouts);

       end
       
       
       function node = getnode(obj,idx)
           valueStruct = struct('individual', obj.getIndividual(idx), ...
                                'knockouts', obj.knockoutsStructure(idx), ...
                                'wins', obj.winsStructure(idx), ...
                                'visits', obj.visitsStructure(idx), ...
                                'childrenRemaining', obj.remainingChildrenStructure(idx));
                            
           node = struct('parent', obj.treeStructure(idx), 'ID', idx, 'value', valueStruct, 'childrenIDs', obj.childrenIDsStructure(idx));
       end
       
       function W = getWeights(obj, idx)
           W = obj.weightStructure(idx);
       end
       
       function Win = getInputs(obj, idx)
           Win = obj.inputStructure(idx);
       end
       
       function Wout = getOutputs(obj, idx)
           Wout = obj.outputStructure(idx);
           Wout = Wout{1};
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
       

       function setWeights(obj, idx, val) 
           obj.weightStructure(idx) = val;
       end
       
       function setInputs(obj, idx, val)
           obj.inputStructure(idx) = val;
       end
       
       function setOutputs(obj, idx, val)
           obj.outputStructure(idx) = {val};
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
       end

       function leaf = isleaf(obj, idx)
           leaf = isempty(obj.getchildren(idx));
       end
       
   end
end

% node format: {struct: parentID | ID | value | [childrenIDs] }
% tree format: {vector: parentID | parentID | parentID ... }


