clear
cost = 0.002;

% Burg entropy
load('../Result/Daily&Monthly_Burg(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);

cnt = 1;
for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;                
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;            
            end
        end
        
        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        cnt = cnt + 1;
    end
end
wealth_daily_cost_burg = wealth_robust_daily_cost;


% x2 distance
load('../Result/Daily&Monthly_x2(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);

cnt = 1;

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;                
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;            
            end
        end
        
        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        cnt = cnt + 1;
    end
end

wealth_daily_cost_x2 = wealth_robust_daily_cost;


% Modified x2 distance
load('../Result/Daily&Monthly_Modified_x2(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);

cnt = 1;

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
               
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;            
            end
        end


        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        cnt = cnt + 1;
    end
end

wealth_daily_cost_modified_x2 = wealth_robust_daily_cost;


% Hellinger distance
load('../Result/Daily&Monthly_Hellinger(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);

cnt = 1;

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;                
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;            
            end
        end
        
        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        cnt = cnt + 1;
    end
end

wealth_daily_cost_hellinger = wealth_robust_daily_cost;


% Variation distance
load('../Result/Daily&Monthly_Variation(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);

cnt = 1;

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;                
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;            
            end
        end
        
        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        cnt = cnt + 1;
    end
end

wealth_daily_cost_variation = wealth_robust_daily_cost;


% KL
load('../Result/Daily&Monthly_KL(d=0.1).mat')
return_robust = zeros(test_length * p, 1);
return_classical = zeros(test_length * p, 1);
return_n = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);
wealth_classical_daily_cost = eye(1, test_length * p + 1);
wealth_n_daily_cost = eye(1, test_length * p + 1);
cnt = 1;

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                cost_classical = norm(portfolio_classical(i, :), 1);                        
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1 - cost_classical * cost;
                return_n(cnt) = wealth_n_daily(cnt + 1) / wealth_n_daily(cnt) - 1;
                
            else          
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1;
                return_n(cnt) = wealth_n_daily(cnt + 1) / wealth_n_daily(cnt) - 1;
            end
        else
            if j == 1            
                cost_robust = norm(portfolio_robust(i, :) - portfolio_robust(i - 1, :), 1);
                cost_classical = norm(portfolio_classical(i, :) - portfolio_classical(i - 1, :), 1);                
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1 - cost_classical * cost;
                return_n(cnt) = wealth_n_daily(cnt + 1) / wealth_n_daily(cnt) - 1;
                
            else           
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1;
                return_n(cnt) = wealth_n_daily(cnt + 1) / wealth_n_daily(cnt) - 1;
            
            end
        end
        
 
        wealth_robust_daily_cost(cnt + 1) = wealth_robust_daily_cost(cnt) * (return_robust(cnt) + 1);
        wealth_classical_daily_cost(cnt + 1) = wealth_classical_daily_cost(cnt) * (return_classical(cnt) + 1);
        wealth_n_daily_cost(cnt + 1) = wealth_n_daily_cost(cnt) * (return_n(cnt) + 1);
        cnt = cnt + 1;
    end
end

wealth_daily_cost_kl = wealth_robust_daily_cost;

relative_daily_cost_burg = wealth_daily_cost_burg ./ wealth_classical_daily_cost - 1;
relative_daily_cost_x2 = wealth_daily_cost_x2 ./ wealth_classical_daily_cost - 1;
relative_daily_cost_modified_x2 = wealth_daily_cost_modified_x2 ./ wealth_classical_daily_cost - 1;
relative_daily_cost_hellinger = wealth_daily_cost_hellinger ./ wealth_classical_daily_cost - 1;
relative_daily_cost_variation = wealth_daily_cost_variation ./ wealth_classical_daily_cost - 1;
relative_daily_cost_kl = wealth_daily_cost_kl ./ wealth_classical_daily_cost - 1;

figure(1)
time_array = datenum(time_new);
plot(time_array, wealth_daily_cost_burg, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_cost_x2, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_cost_modified_x2, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_cost_hellinger, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_cost_variation, 'LineWidth', 1);
hold on
plot(time_array, wealth_daily_cost_kl, 'LineWidth', 1);
hold on
plot(time_array, wealth_classical_daily_cost, 'LineWidth', 1);
legend('Burg entropy', 'x2', 'Modified x2', 'Hellinger', 'Variation', 'KL', 'Classical', 'Location', 'northwest');
xlabel('Time');
ylabel('Wealth');
datetick('x', 25);
xlim([time_array(1), time_array(481)]);
title('Daily Wealth After Transaction Cost');


figure(2)
time_array = datenum(time_new); 
plot(time_array, relative_daily_cost_burg, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_cost_x2, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_cost_modified_x2, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_cost_hellinger, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_cost_variation, 'LineWidth', 1);
hold on
plot(time_array, relative_daily_cost_kl, 'LineWidth', 1);

legend('Burg entropy', 'x2', 'Modified x2', 'Hellinger', 'Variation', 'KL', 'Location', 'southwest');
xlabel('Time');
ylabel('Relative Wealth Ratio')
datetick('x', 25);
xlim([time_array(1), time_array(481)]);
title('Relative Daily Wealth of Robust/Classical Portfolios After Transaction Cost');



