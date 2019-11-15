clear
load('../Result/Daily&Monthly_Burg(d=0.1).mat')
wealth_daily_burg = wealth_robust_daily;

load('../Result/Daily&Monthly_x2(d=0.1).mat')
wealth_daily_x2 = wealth_robust_daily;

load('../Result/Daily&Monthly_Modified_x2(d=0.1).mat')
wealth_daily_modified_x2 = wealth_robust_daily;

load('../Result/Daily&Monthly_Hellinger(d=0.1).mat')
wealth_daily_hellinger = wealth_robust_daily;

load('../Result/Daily&Monthly_Variation(d=0.1).mat')
wealth_daily_variation = wealth_robust_daily;

load('../Result/Daily&Monthly_KL(d=0.1).mat')
wealth_daily_kl = wealth_robust_daily;

relative_daily_burg = wealth_daily_burg ./ wealth_classical_daily - 1;
relative_daily_x2 = wealth_daily_x2 ./ wealth_classical_daily - 1;
relative_daily_modified_x2 = wealth_daily_modified_x2 ./ wealth_classical_daily - 1;
relative_daily_hellinger = wealth_daily_hellinger ./ wealth_classical_daily - 1;
relative_daily_variation = wealth_daily_variation ./ wealth_classical_daily - 1;
relative_daily_kl = wealth_daily_kl ./ wealth_classical_daily - 1;

figure(1)
time_array = datenum(time_new);
plot(time_array, wealth_daily_burg, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_x2, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_modified_x2, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_hellinger, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_variation, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_kl, 'LineWidth', 1);
hold on
plot(time_array, wealth_classical_daily, 'LineWidth', 1);
legend('Burg entropy', 'x2', 'Modified x2', 'Hellinger', 'Variation', 'KL', 'Classical', 'Location', 'northwest');
xlabel('Time');
ylabel('Wealth');
datetick('x', 25);
xlim([time_array(1), time_array(end)]);
title('Daily Wealth Under Several Strategies');


figure(2)
time_array = datenum(time_new); 
plot(time_array, relative_daily_burg, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_x2, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_modified_x2, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_hellinger, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_variation, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_kl, 'LineWidth', 1);

legend('Burg entropy', 'x2', 'Modified x2', 'Hellinger', 'Variation', 'KL', 'Location', 'southwest');
xlabel('Time');
ylabel('Relative Wealth Ratio')
datetick('x', 25);
xlim([time_array(1), time_array(end)]);
title('Relative Daily Wealth of Robust/Classical Portfolios');
