clear
%load('../Result/Real_KL(d=0.1).mat')
%load('../Result/Real_Burg(d=0.1).mat')
%load('../Result/Real_x2(d=0.1).mat')
load('../Result/Real_Modified_x2(d=0.1).mat')
%load('../Result/Real_Hellinger(d=0.1).mat')
%load('../Result/Real_Variation(d=0.1).mat')

[data_return, time_cell] = xlsread('../Data/data.xlsx');
[N, n] = size(data_return);
p = 20;
train_length = 45;
test_length = 24;
time = string(time_cell((train_length) * p + 1 : (test_length + train_length) * p + 1));
time_new = datetime(time,'InputFormat','yyyy-MM-dd');
T = p * train_length;

wealth_robust_daily = eye(1, test_length * p + 1);
wealth_classical_daily = eye(1, test_length * p + 1);
wealth_n_daily = eye(1, test_length * p + 1);

portfolio_robust;
x_n = ones(n, 1) / n;
i = 1;

for period = 1 : test_length
    return_test = data_return((period + train_length - 1) * p + 1 : (period + train_length) * p, :);
    
    for day = 1 : p
        return_daily = return_test(day, :);
        return_robust = (1 + return_daily) * portfolio_robust(period, :)';
        return_classical = (1 + return_daily) * portfolio_classical(period, :)';
        return_n = (1 + return_daily) * x_n;
        
        wealth_robust_daily(i + 1) = wealth_robust_daily(i) * return_robust;
        wealth_classical_daily(i + 1) = wealth_classical_daily(i) * return_classical;
        wealth_n_daily(i + 1) = wealth_n_daily(i) * return_n;
        
        i = i + 1;
    end
    
end

figure
subplot(2, 1, 1)
time_array = datenum(time_new);
txt = ['Daily Wealth Under Three Strategies (d = ', num2str(d), ')'];
plot(time_array, wealth_robust_daily);
hold on
plot(time_array, wealth_classical_daily);
hold on
plot(time_array, wealth_n_daily);
legend('robust', 'classical', '1/n', 'Location', 'northwest');
xlabel('Time');
ylabel('Wealth')
datetick('x', 25)
xlim([time_array(1), time_array(481)]);
title(txt);   


subplot(2, 1, 2)
txt = ['Monthly Wealth Under Three Strategies (d = ', num2str(d), ')'];
plot(t_array, wealth_robust);
hold on
plot(t_array, wealth_classical);
hold on
plot(t_array, wealth_n);
legend('robust', 'classical', '1/n', 'Location', 'northwest');
xlabel('Period');
ylabel('Wealth');
title(txt);      
