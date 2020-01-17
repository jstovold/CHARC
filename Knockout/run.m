
if false
    clear
    close all

    [all_diffs_25, metrics_a_25, metrics_b_25, ks_outs_25, atest_outs_25, pop_25] = run_full_tests(250,25);
    [all_diffs_50, metrics_a_50, metrics_b_50, ks_outs_50, atest_outs_50, pop_50] = run_full_tests(250,50);
    [all_diffs_100, metrics_a_100, metrics_b_100, ks_outs_100, atest_outs_100, pop_100] = run_full_tests(250,100);
    [all_diffs_200,  metrics_a_200, metrics_b_200, ks_outs_200, atest_outs_200, pop_200] = run_full_tests(250,200);

end

scatter(metrics_a_25(1,:), metrics_a_25(2,:))



subplot(2,3,1)
scatter(metrics_a_25(1,:), metrics_a_25(2,:), 'rx')
hold on
scatter(metrics_b_25(1,:), metrics_b_25(2,:), 'bo')
xlabel('KR');
ylabel('GR');
hold off

subplot(2,3,2)
scatter(metrics_a_25(2,:), metrics_a_25(3,:), 'rx')
hold on
scatter(metrics_b_25(2,:), metrics_b_25(3,:), 'bo')
xlabel('GR');
ylabel('MC');
hold off

subplot(2,3,3)
scatter(metrics_a_25(1,:), metrics_a_25(3,:), 'rx')
hold on
scatter(metrics_b_25(1,:), metrics_b_25(3,:), 'bo')
xlabel('KR')
ylabel('MC')
hold off





subplot(2,3,4)
scatter(metrics_a_50(1,:), metrics_a_50(2,:), 'rx')
hold on
scatter(metrics_b_50(1,:), metrics_b_50(2,:), 'bo')
xlabel('KR');
ylabel('GR');
hold off

subplot(2,3,5)
scatter(metrics_a_50(2,:), metrics_a_50(3,:), 'rx')
hold on
scatter(metrics_b_50(2,:), metrics_b_50(3,:), 'bo')
xlabel('GR');
ylabel('MC');
hold off

subplot(2,3,6)
scatter(metrics_a_50(1,:), metrics_a_50(3,:), 'rx')
hold on
scatter(metrics_b_50(1,:), metrics_b_50(3,:), 'bo')
xlabel('KR')
ylabel('MC')
hold off








subplot(1,3,1)
scatter(metrics_a_50(1,:), metrics_a_50(2,:), 'bo')
hold on
scatter(metrics_b_50(1,:), metrics_b_50(2,:), 'rx')
hold off

subplot(1,3,2)
scatter(metrics_a_50(2,:), metrics_a_50(3,:), 'bo')
hold on
scatter(metrics_b_50(2,:), metrics_b_50(3,:), 'rx')
hold off

subplot(1,3,3)
scatter(metrics_a_50(1,:), metrics_a_50(3,:), 'bo')
hold on
scatter(metrics_b_50(1,:), metrics_b_50(3,:), 'rx')
hold off




subplot(1,3,1)
scatter(metrics_a_100(1,:), metrics_a_100(2,:), 'bo')
hold on
scatter(metrics_b_100(1,:), metrics_b_100(2,:), 'rx')
hold off

subplot(1,3,2)
scatter(metrics_a_100(2,:), metrics_a_100(3,:), 'bo')
hold on
scatter(metrics_b_100(2,:), metrics_b_100(3,:), 'rx')
hold off

subplot(1,3,3)
scatter(metrics_a_100(1,:), metrics_a_100(3,:), 'bo')
hold on
scatter(metrics_b_100(1,:), metrics_b_100(3,:), 'rx')
hold off


subplot(1,3,1)
scatter(metrics_a_200(1,:), metrics_a_200(2,:), 'bo')
hold on
scatter(metrics_b_200(1,:), metrics_b_200(2,:), 'rx')
hold off

subplot(1,3,2)
scatter(metrics_a_200(2,:), metrics_a_200(3,:), 'bo')
hold on
scatter(metrics_b_200(2,:), metrics_b_200(3,:), 'rx')
hold off

subplot(1,3,3)
scatter(metrics_a_200(1,:), metrics_a_200(3,:), 'bo')
hold on
scatter(metrics_b_200(1,:), metrics_b_200(3,:), 'rx')
hold off

