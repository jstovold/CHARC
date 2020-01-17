%%
function [all_diffs, all_metrics_a, all_metrics_b, ks_outs, atest_outs, population] = run_full_tests(pop_size,res_size)
    
    [population, config] = create_population(pop_size, res_size);
    close all
    
    all_diffs = zeros(4,pop_size);
    all_metrics_a = zeros(4,pop_size);
    all_metrics_b = zeros(4,pop_size);
    for p = 1:length(population)
        [diff, metrics_a, metrics_b] = run_knockout(population(p), config);
        all_diffs(:,p) = diff;
        all_metrics_a(:,p) = metrics_a;
        all_metrics_b(:,p) = metrics_b;
        disp(p);
    end
    
    ks_outs = zeros(4,2);
    atest_outs = zeros(4,1);
    [h,p,~] = kstest2(all_metrics_a(1,:)', all_metrics_b(1,:)', 'Tail', 'Larger');
    ks_outs(1,:) = [h,p];
    a = Atest(all_metrics_a(1,:)', all_metrics_b(1,:)');
    atest_outs(1) = a;
    
    [h,p,~] = kstest2(all_metrics_a(2,:)', all_metrics_b(2,:)', 'Tail', 'Smaller');
    ks_outs(2,:) = [h,p];
    a = Atest(all_metrics_a(2,:)', all_metrics_b(2,:)');
    atest_outs(2) = a;
    
    [h,p,~] = kstest2(all_metrics_a(3,:)', all_metrics_b(3,:)', 'Tail', 'Larger');
    ks_outs(3,:) = [h,p];
    a = Atest(all_metrics_a(3,:)', all_metrics_b(3,:)');
    atest_outs(3) = a;
    
    [h,p,~] = kstest2(all_metrics_a(4,:)', all_metrics_b(4,:)', 'Tail', 'Smaller');
    ks_outs(4,:) = [h,p];
    a = Atest(all_metrics_a(4,:)', all_metrics_b(4,:)');
    atest_outs(4) = a;
    
    combined = [all_metrics_a(1,:); all_metrics_b(1,:); all_metrics_a(2,:); all_metrics_b(2,:); all_metrics_a(3,:); all_metrics_b(3,:); all_metrics_a(4,:); all_metrics_b(4,:)];
    positions = [1 1.25 2 2.25 3 3.25 4 4.25];
    boxplot(combined', 'positions', positions);
    
    boxplot(all_diffs')
    
    
    