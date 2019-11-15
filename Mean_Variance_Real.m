clear
cnt = 1;
cost = 0.002;
%load('../Result/Daily&Monthly_KL(d=0.1).mat')
%load('../Result/Daily&Monthly_Burg(d=0.1).mat')
%load('../Result/Daily&Monthly_x2(d=0.1).mat')
load('../Result/Daily&Monthly_Modified_x2(d=0.1).mat')
%load('../Result/Daily&Monthly_Hellinger(d=0.1).mat')
%load('../Result/Daily&Monthly_Variation(d=0.1).mat')

return_robust = zeros(test_length * p, 1);
return_classical = zeros(test_length * p, 1);
return_n = zeros(test_length * p, 1);
wealth_robust_daily_cost = eye(1, test_length * p + 1);
wealth_classical_daily_cost = eye(1, test_length * p + 1);
wealth_n_daily_cost = eye(1, test_length * p + 1);

for i = 1 : test_length
    for j = 1 : p
        if i == 1 
            if j == 1
                cost_robust = norm(portfolio_robust(i, :), 1);
                cost_classical = norm(portfolio_classical(i, :), 1);                        
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost * 10 / 9;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1 - cost_classical * cost* 10 / 9;
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
                return_robust(cnt) = wealth_robust_daily(cnt + 1) / wealth_robust_daily(cnt) - 1 - cost_robust * cost * 10 / 9;
                return_classical(cnt) = wealth_classical_daily(cnt + 1) / wealth_classical_daily(cnt) - 1 - cost_classical * cost * 10 / 9;
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

% return, volatility, sharpe ratio
annual_return_robust = (wealth_robust_daily_cost(end) / wealth_robust_daily_cost(1) - 1) / (test_length * p) * 250;
volatility_robust = std(return_robust) * sqrt(250);
sharpe_robust = annual_return_robust / volatility_robust;

annual_return_classical = (wealth_classical_daily_cost(end) / wealth_classical_daily_cost(1) - 1) / (test_length * p) * 250;
volatility_classical = std(return_classical) * sqrt(250);
sharpe_classical = annual_return_classical / volatility_classical;

annual_return_n = (wealth_n_daily_cost(end) / wealth_n_daily_cost(1) - 1) / (test_length * p) * 250;
volatility_n = std(return_n) * sqrt(250);
sharpe_n = annual_return_n / volatility_n;

% beta, alpha
cov_matrix_robust = cov(return_robust, return_n);
cov_robust = cov_matrix_robust(1, 2);
beta_robust = cov_robust / var(return_n);
alpha_robust = annual_return_robust - beta_robust * annual_return_n;

cov_matrix_classical = cov(return_classical, return_n);
cov_classical = cov_matrix_classical(1, 2);
beta_classical = cov_classical / var(return_n);
alpha_classical = annual_return_classical - beta_classical * annual_return_n;

% max drawdown
drawdown_robust = zeros(test_length * p, 1);
drawdown_classical = zeros(test_length * p, 1);
drawdown_n = zeros(test_length * p, 1);

for i = 1 : test_length * p
    drawdown_robust(i) = 1 - wealth_robust_daily_cost(i) / max(wealth_robust_daily_cost(1 : i));
    drawdown_classical(i) = 1 - wealth_classical_daily_cost(i) / max(wealth_classical_daily_cost(1 : i));
    drawdown_n(i) = 1 - wealth_n_daily_cost(i) / max(wealth_n_daily_cost(1 : i));
end

max_drawdown_robust = max(drawdown_robust);
max_drawdown_classical = max(drawdown_classical);
max_drawdown_n = max(drawdown_n);