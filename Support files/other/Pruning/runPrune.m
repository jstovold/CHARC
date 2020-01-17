%%
function pruned_database = runPrune(database, config)

if isempty(gcp) && config.parallel
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

rng(1,'twister')
num_iter = 4;
ppm = ParforProgMon('Database complete: ', num_iter);
for indx = 1:num_iter
    warning('off','all')
    individual = database(indx);
    for p = 1:5
        %[individual,individual.behaviours,~,error] = maxPruning(@getMetrics, 'behaviours', individual, database(indx).behaviours, [0 0 0.5], config);
        [individual,individual.behaviours,~,error] = maxPruning(@getMetrics, 'behaviours', individual, getMetrics(database(indx), config), [0 0 0.5], config);
        fprintf('Node: %d, iter: %d, error: %d \n',indx,p,full(error))
    end
    
    pruned_database(indx) = individual;
%     ppm.increment();
end

% figure1 = figure;

% for plot_indx = 1:length(pruned_database)

% M = pruned_database(plot_indx).W{1,1};
% [MATreordered_b] = reorder_matrix(M~=0,'line',0);
% subplot(1,2,1)
% imagesc(M)
% colormap(gca,bluewhitered)
% title(num2str(pruned_database(plot_indx).behaviours))
% 
% [MATreordered_W] = reorder_matrix(abs(M),'line',0);
% subplot(1,2,2)
% imagesc(MATreordered_W)
% colormap(gca,bluewhitered)
% title(num2str(nnz(M)))
% drawnow
% pause(0.1)

% end


