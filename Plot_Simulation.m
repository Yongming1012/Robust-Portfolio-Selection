load('../Result/Modified_x2_SR.mat')
sharpe_mean_Modified_x2 = sharpe_array_robust_classical;
sharpe_worst_Modified_x2 = sharpe_array_worst_robust_classical;

load('../Result/Burg_SR.mat')
sharpe_mean_Burg = sharpe_array_robust_classical;
sharpe_worst_Burg = sharpe_array_worst_robust_classical;

load('../Result/KL_SR.mat')
sharpe_mean_KL = sharpe_array_robust_classical;
sharpe_worst_KL = sharpe_array_worst_robust_classical;

load('../Result/Hellinger_SR.mat')
sharpe_mean_Hellinger = sharpe_array_robust_classical;
sharpe_worst_Hellinger = sharpe_array_worst_robust_classical;


d_array = [0: 0.02: 1];
x = [0: 0.01: 1];
sharpe_mean_Modified_x2_y = spline(d_array, sharpe_mean_Modified_x2, x);
sharpe_worst_Modified_x2_y = spline(d_array, sharpe_worst_Modified_x2, x);
sharpe_mean_Burg_y = spline(d_array, sharpe_mean_Burg, x);
sharpe_worst_Burg_y = spline(d_array, sharpe_worst_Burg, x);
sharpe_mean_KL_y = spline(d_array, sharpe_mean_KL, x);
sharpe_worst_KL_y = spline(d_array, sharpe_worst_KL, x);
sharpe_mean_Hellinger_y = spline(d_array, sharpe_mean_Hellinger, x);
sharpe_worst_Hellinger_y = spline(d_array, sharpe_worst_Hellinger, x);

figure(1)
plot(d_array, sharpe_mean_Modified_x2, 'LineWidth', 2);
hold on
plot(d_array, sharpe_mean_Burg, 'LineWidth', 2);
hold on
plot(d_array, sharpe_mean_KL, 'LineWidth', 2);
hold on
plot(d_array, sharpe_mean_Hellinger, 'LineWidth', 2);
legend('Modified x2 distance', 'Burg entropy', 'KL', 'Hellinger distance');
xlabel('d');
ylabel('Robust/Classical Sharpe Ratio')
title('Robust/Classical Mean Sharpe Ratio');

figure(2)
plot(d_array, sharpe_worst_Modified_x2, 'LineWidth', 2);
hold on
plot(d_array, sharpe_worst_Burg, 'LineWidth', 2);
hold on
plot(d_array, sharpe_worst_KL, 'LineWidth', 2);
hold on
plot(d_array, sharpe_worst_Hellinger, 'LineWidth', 2);
legend('Modified_x2_distance', 'Burg entropy', 'KL', 'Hellinger distance');
xlabel('d');
ylabel('Robust/Classical Sharpe Ratio')
title('Robust/Classical Worst Sharpe Ratio');