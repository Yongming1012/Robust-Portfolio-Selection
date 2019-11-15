clear all
warning('off');

data_return = xlsread('../Data/data.xlsx');
eta = xlsread('../Data/eta_real.xlsx');


% n: number of assets  n_sample: time length
[n_sample, n] = size(data_return); 
p = 20;
a = 0.001; % minimal return
d = 0.1;

train_length = 45;
test_length = 24;
T = p * train_length;

wealth_robust = eye(1, test_length + 1);
wealth_classical = eye(1, test_length + 1);
wealth_n = eye(1, test_length + 1);

relative_cost = zeros(1, test_length - 1);
portfolio_robust = zeros(test_length, n);
portfolio_classical = zeros(test_length, n);

return_ori = 0;
step_size_array = [0.01: 0.04: 0.4];
mu = 0.5;

for period = 1 : test_length
    
    R = data_return((period - 1) * p + 1 : (period + train_length - 1) * p, :)';
    return_test = data_return((period + train_length - 1) * p + 1 : (period + train_length) * p, :)';
    
    %-------------------------------------------------------%
    % classical strategy
    disp('Begin to solve the classical portfolio problem...');
    
    % Classical maximum variance problem
    cov_classical = zeros(n, n);
    for k = 1 : T
        r = R(:, k);
        cov_classical = cov_classical + eta(k) * (r - R * eta) * (r - R * eta)';
    end
    cov = (cov_classical + cov_classical') / 2;
    cov_classical = cov;
    
    %{
    cvx_begin quiet
    cvx_solver mosek
    variables x_classical(n, 1)
    minimize( 1 )
    subject to
        (R * eta)' * x_classical >= a;
        ones(1, n) * x_classical == 1;
        x_classical >= 0;
    cvx_end

    tolerance = 1e-14;
    error = 1;
    obj_array = [10];
    count = 0;
    for iteration = 1 : 1000
        iteration
        gradient = 2 * cov_classical' * x_classical;
        step_size = 0.1;
        x_poj_ori = x_classical - step_size * gradient;

        cvx_begin quiet
        cvx_solver mosek
        variables x_poj(n, 1)
        minimize( norm(x_poj - x_poj_ori, 2) )
        subject to
            (R * eta)' * x_poj >= a;
            ones(1, n) * x_poj == 1;
            x_poj >= 0;
        cvx_end
        obj = x_poj' * cov_classical * x_poj;
        error = abs(obj_array(iteration) - obj);
        error_iteration = obj_array(iteration) - obj;

        if error < tolerance
            break
        end
        if error_iteration <= 0
            count = count + 1;
        end
        if count == 2
            break
        end

        x_classical = x_poj;
        obj_array = [obj_array, obj];
        error
        error_iteration
        obj_array(iteration + 1)

    end
    %}
    cvx_begin quiet
    cvx_solver mosek
    variables x_classical(n, 1)
    minimize( x_classical' * cov_classical * x_classical )
    subject to
        (R * eta)' * x_classical >= a;
        ones(1, n) * x_classical == 1;
        x_classical >= 0;
    cvx_end
    
    var_classical = x_classical' * cov_classical * x_classical
    portfolio_classical(period, :) = x_classical;
    
    %--------------------------------------------%
    % robust strategy
    disp('Begin to solve the robust portfolio problem...');
    
    % Initialization: find the feasible solution
    cvx_begin quiet
    cvx_solver mosek
    variables x(n, 1) alpha_rps gamma_rps y(T, 1) 
    minimize( 1 )
    ones(1, n) * x == 1;
    eta' * y + alpha_rps * d + gamma_rps <= -a;
    max(-y' * eta, -alpha_rps) <= (R * eta)' * x + gamma_rps;
    alpha_rps * ones(T, 1) + y >= 0;
    x >= 0;
    gamma_rps >= 0;
    alpha_rps >= 0;
    cvx_end
    
    if d == 0
        pi = eta;
    else
        % find the optimal pi
        w = R' * x;
        cvx_begin quiet
        cvx_solver mosek
        variables pi_opt(T, 1)
        maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
        sum(abs(pi_opt - eta)) <= d;
        ones(1, T) * pi_opt == 1;
        pi_opt >= 0;
        cvx_end

        pi = pi_opt;
        for k = 1 : T
            if pi(k) < 0
                pi(k) = 0.0;
            end
        end
    end

    covariance = zeros(n, n);
    for k = 1 : T
        r = R(:, k);
        covariance = covariance + pi(k) * (r - R * pi) * (r - R * pi)';
    end
    cov = (covariance + covariance') / 2;
    covariance = cov;
    %{
    cvx_begin quiet
    cvx_solver mosek
    variables x(n, 1) alpha_rps gamma_rps y(T, 1) 
    minimize( x' * covariance * x )
    ones(1, n) * x == 1;
    eta' * y + alpha_rps * d + gamma_rps <= -a;
    max(-y' * eta, -alpha_rps) <= (R * eta)' * x + gamma_rps;
    alpha_rps * ones(T, 1) + y >= 0;
    x >= 0;
    gamma_rps >= 0;
    alpha_rps >= 0;
    cvx_end
    
    % find the optimal pi
    w = R' * x;
    cvx_begin quiet
    cvx_solver mosek
    variables pi_opt(T, 1)
    maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
    sum(abs(pi_opt - eta)) <= d;
    ones(1, T) * pi_opt == 1;
    pi_opt >= 0;
    cvx_end

    pi_new = pi_opt;
    for k = 1 : T
        if pi_new(k) < 0
            pi_new(k) = 0.0;
        end
    end
    
    covariance = zeros(n, n);
    for k = 1 : T
        r = R(:, k);
        covariance = covariance + pi_new(k) * (r - R * pi_new) * (r - R * pi_new)';
    end
    cov = (covariance + covariance') / 2;
    covariance = cov;
    %}
   
    x' * covariance * x

    obj_array = [x' * covariance * x];
    error_array = [10];

    tolerance = 1e-9;
    error = 1;
    count = 0;

    for iteration = 1 : 500
        iteration

        gradient = 2 * covariance' * x;
        
        
        for step = 1 : 10
            step_size_tmp = step_size_array(11 - step);
            x_tmp = x - step_size_tmp * gradient;
            
            cvx_begin quiet
            cvx_solver mosek
            variables x_poj(n, 1) alpha_rps gamma_rps y(T, 1) 
            minimize( norm(x_poj - x_tmp, 2) )
            ones(1, n) * x_poj == 1;
            eta' * y + alpha_rps * d + gamma_rps <= -a;
            max(-y' * eta, -alpha_rps) <= (R * eta)' * x_poj + gamma_rps;
            alpha_rps * ones(T, 1) + y >= 0;
            x_poj >= 0;
            gamma_rps >= 0;
            alpha_rps >= 0;
            cvx_end
            
            if any(isnan(x_poj)) > 0
                continue
            end
            
            f_x_step = x_poj' * covariance * x_poj;
            f_x = x' * covariance * x;
            
            if (f_x - f_x_step) >= mu * gradient' * (x - x_poj)
                step_size = step_size_tmp;
                x_new = x_poj;
                break
            else
                continue
            end

        end
        
        
        step_size
        
        
        %step_size = 0.01;
        %x_poj_ori = x - step_size * gradient;
        
        %{
        cvx_begin quiet
        cvx_solver mosek
        variables x_poj(n, 1) alpha_rps gamma_rps y(T, 1) 
        minimize( norm(x_poj - x_poj_ori, 2) )
        ones(1, n) * x_poj == 1;
        eta' * y + alpha_rps * d + gamma_rps <= -a;
        max(-y' * eta, -alpha_rps) <= (R * eta)' * x_poj + gamma_rps;
        alpha_rps * ones(T, 1) + y >= 0;
        x_poj >= 0;
        gamma_rps >= 0;
        alpha_rps >= 0;
        cvx_end
        %}
        
        %direction = x_poj - x;
        
        %{
        c = 0.7;
        cvx_begin quiet
        cvx_solver mosek
        variables lambda_k
        maximize( lambda_k )
        f_k_lambda = quad_form(x + lambda_k * direction, covariance);
        f_k = quad_form(x, covariance);
        f_k_lambda - f_k < c * lambda_k * gradient' * direction;
        lambda_k > 0;
        lambda_k < 1;
        cvx_end
        %}
        
        
        %f_k_lambda - f_k
        
        %step_size = lambda_k
        %step_size = 1 / iteration;

        %x_new = x_poj;
        
        if d == 0
            pi_new = eta;
        else
            % find the optimal pi
            w = R' * x_new;
            cvx_begin quiet
            cvx_solver mosek
            variables pi_opt(T, 1)
            maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
            sum(abs(pi_opt - eta)) <= d;
            ones(1, T) * pi_opt == 1;
            pi_opt >= 0;
            cvx_end
            if any(isnan(pi_opt)) > 0
                pi_new = pi;
            else
                pi_new = pi_opt;
            end
            for k = 1 : T
                if pi_new(k) < 0
                    pi_new(k) = 0.0;
                end
            end        
        end
        
        cov_new = zeros(n, n);
        for k = 1 : T
            r = R(:, k);
            cov_new = cov_new + pi_new(k) * (r - R * pi_new) * (r - R * pi_new)';
        end
        cov = (cov_new + cov_new') / 2;
        cov_new = cov;

        obj = x_new' * cov_new * x_new;
        error = abs(obj_array(iteration) - obj);
        %step = norm(x_new - x);        
        error_iteration = obj_array(iteration) - obj

        if floor(log10(x_new' * cov_classical * x_new)) < floor(log10(x_classical' * cov_classical * x_classical))
            disp('Variance exceeds!');
            break
        end
        
        if error <= tolerance || (x_new' * cov_classical * x_new) < 0 || obj < 0
            break
        else
            x = x_new;
            pi = pi_new;
            covariance = cov_new;
        end
        
        if error_iteration <= 0
            count = count + 1;
        end
        if count == 2
            break
        end
        
        error_array = [error_array, error_iteration];
        %if error_array(iteration) < error_array(iteration + 1)
        %    break
        %end    
        obj_array = [obj_array, obj];
        fprintf('Current worst return: %e \n', (R * pi)' * x);
        fprintf('Current mean return: %e \n', (R * eta)' * x); 
        fprintf('Current worst variance: %e \n', obj_array(iteration + 1));
        fprintf('Current mean variance: %e \n', x' * cov_classical * x);

    end
    disp('Robust portfolio problem solved.');
    portfolio_robust(period, :) = x;
    %--------------------------------------------------------------%

    % 1 / n strategy
    x_n = ones(n, 1) / n;
    
    %--------------------------------------------------------------%
    % calculate relative transaction cost
    if period ~= 1
        cost_robust = norm(portfolio_robust(period, :) - portfolio_robust(period - 1, :), 1);
        cost_classical = norm(portfolio_classical(period, :) - portfolio_classical(period - 1, :), 1);
        relative_cost(period - 1) = cost_robust / cost_classical;
    end
    
    %--------------------------------------------------------------%
    % calculate return
    total_return = ones(n, 1);
    for j = 1 : p
        return_today = return_test(:, j);
        total_return = (1 + return_today) .* total_return;
    end
    
    wealth_robust(period + 1) = (total_return' * x) * wealth_robust(period);
    wealth_classical(period + 1) = (total_return' * x_classical) * wealth_classical(period);
    wealth_n(period + 1) = (total_return' * x_n) * wealth_n(period);
    
    fprintf('Period: %d \n', period)
    fprintf('My strategy: %f \n', wealth_robust(period + 1))
    fprintf('Classical statrgy: %f \n', wealth_classical(period + 1))
    fprintf('1/n strategy: %f \n', wealth_n(period + 1))
    fprintf('\n')
end

t_array = [0 : 1 : test_length];
cost_array = [2 : 1 : test_length];


figure(1)
txt = ['Wealth Under Three Strategies (d = ', num2str(d), ')'];
plot(t_array, wealth_robust);
hold on
plot(t_array, wealth_classical);
hold on
plot(t_array, wealth_n);
legend('robust', 'classical', '1/n');
xlabel('Period');
ylabel('Wealth');
title(txt);      


figure(2)
average_cost = mean(relative_cost);
txt_1 = ['Average Cost = ', num2str(average_cost)];
txt_2 = ['Relative Cost of Robust vs. Classical Portfolios (d = ', num2str(d), ')'];
plot(cost_array, relative_cost);
xlabel('Period');
ylabel('Relative Cost');
text(16, 13, txt_1);
title(txt_2);

disp(['Mean return of robust portfolio: ', num2str(mean(wealth_robust))]);
disp(['Variance of robust portfolio: ', num2str(var(wealth_robust))]);
disp(['Mean return of classical portfolio: ', num2str(mean(wealth_classical))]);
disp(['Variance return of classical portfolio: ', num2str(var(wealth_classical))]);
disp(['Mean return of 1/n performance: ', num2str(mean(wealth_n))]);
disp(['Variance of 1/n performance: ', num2str(var(wealth_n))]);

save('../Result/Real_Variation(d=0.1).mat', 'd', 'average_cost', 'cost_array', 'portfolio_classical', 'portfolio_robust', 'relative_cost', 't_array', 'wealth_classical', 'wealth_n', 'wealth_robust');
