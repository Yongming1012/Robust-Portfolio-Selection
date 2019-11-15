clear all
warning('off');

n = 500; % number of assets
T = 50; % time length
a = 2; % minimal return
r_f = 2;

eta_ori = rand(T, 1);
eta_sum = sum(eta_ori);
eta = eta_ori / eta_sum;
eta_2 = eta .* eta;

R = zeros(n, T);
for k = 1 : T
    r = randn(n, 1) + rand(n, 1) * 6 - 1;
    R(:, k) = r;
end

d_array = [0: 0.02: 1];
[n_1, n_d] = size(d_array);
sharpe_array_robust_classical = ones(1, n_d);
sharpe_array_worst_robust_classical = ones(1, n_d);

% Classical maximum variance problem
cov_classical = zeros(n, n);
for k = 1 : T
    r = R(:, k);
    cov_classical = cov_classical + eta(k) * (r - R * eta) * (r - R * eta)';
end
cov = (cov_classical + cov_classical') / 2;
cov_classical = cov;

%-------------------------------------------------------%
% solve the classical problem
disp('Begin to solve the classical portfolio problem...');

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
for iteration = 1 : 500
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
x_classical' * cov_classical * x_classical
%--------------------------------------------%

%{
cvx_begin quiet
cvx_solver mosek
variables x_classical(n, 1)
minimize( x_classical' * cov_classical * x_classical )
subject to
    (R * eta)' * x_classical >= a;
    ones(1, n) * x_classical == 1;
    x_classical >= 0;
cvx_end
%}
sharpe_ratio_classical = ((R * eta - r_f * ones(n, 1))' * x_classical) / sqrt(x_classical' * cov_classical * x_classical);
x_classical' * cov_classical * x_classical
sharpe_ratio_classical
step_size_array = [0.01: 0.02: 0.2];
mu = 0.5;  % c

for d_curr = 1 : n_d
    d = d_array(d_curr);
    d

    % Initialization: find the feasible solution

    cvx_begin quiet
    cvx_solver mosek
    variables x(n, 1) alpha_rps gamma_rps z(T, 1) 
    minimize( 1 )
    ones(1, n) * x == 1;
    gamma_rps + alpha_rps * d + eta' * z <= -a;
    for k = 1 : T
        r = R(:, k);
        r' * x >= rel_entr(alpha_rps, alpha_rps + z(k)) - gamma_rps;
        alpha_rps + z(k) > 0;
    end
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
        sum(kl_div(pi_opt, eta)) <= d;
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
    variables x(n, 1) alpha_rps gamma_rps z(T, 1) 
    minimize( x' * covariance * x )
    ones(1, n) * x == 1;
    gamma_rps + alpha_rps * d + eta' * z <= -a;
    for k = 1 : T
        r = R(:, k);
        r' * x >= rel_entr(alpha_rps, alpha_rps + z(k)) - gamma_rps;
        alpha_rps + z(k) > 0;
    end
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
    sum(kl_div(pi_opt, eta)) <= d;
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
    
    %----------------------------------------------%
   
    obj_array = [x' * covariance * x];
    error_array = [10];

    tolerance = 1e-16;
    error = 1;
    count = 0;

    for iteration = 1 : 1000
        iteration
        
        %{
        % find the optimal pi
        w = R' * x;
        cvx_begin quiet
        cvx_solver mosek
        variables pi_opt(T, 1)
        maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
        sum(kl_div(pi_opt, eta)) <= d;
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
        %}
        
        gradient = 2 * covariance' * x;
        
        
        for step = 1 : 10
            step_size_tmp = step_size_array(11 - step);
            x_tmp = x - step_size_tmp * gradient;
            
            cvx_begin quiet
            cvx_solver mosek
            variables x_poj(n, 1) alpha_rps gamma_rps z(T, 1) 
            minimize( norm(x_poj - x_tmp, 2) )
            ones(1, n) * x_poj == 1;
            gamma_rps + alpha_rps * d + eta' * z <= -a;
            for k = 1 : T
                r = R(:, k);
                r' * x_poj >= rel_entr(alpha_rps, alpha_rps + z(k)) - gamma_rps;
                alpha_rps + z(k) > 0;
            end
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
        variables x_poj(n, 1) alpha_rps gamma_rps z(T, 1) 
        minimize( norm(x_poj - x_poj_ori, 2) )
        ones(1, n) * x_poj == 1;
        gamma_rps + alpha_rps * d + eta' * z <= -a;
        for k = 1 : T
            r = R(:, k);
            r' * x_poj >= rel_entr(alpha_rps, alpha_rps + z(k)) - gamma_rps;
            alpha_rps + z(k) > 0;
        end
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

        x_new = x_poj;
        
        if d == 0
            pi_new = eta;
        else
            % find the optimal pi
            w = R' * x_new;
            cvx_begin quiet
            cvx_solver mosek
            variables pi_opt(T, 1)
            maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
            sum(kl_div(pi_opt, eta)) <= d;
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
            sharpe_ratio_worst_robust = ((R * pi - r_f * ones(n, 1))' * x) / sqrt(obj);
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
        fprintf('Current worst variance: %e \n', obj_array(iteration + 1));
        fprintf('Current mean variance: %e \n', x' * cov_classical * x);
       
    end
    
    (R * eta)' * x
    (R * eta)' * x_classical
    x' * cov_classical * x
    x_classical' * cov_classical * x_classical
    
    sharpe_ratio_robust = real(((R * eta - r_f * ones(n, 1))' * x) / sqrt(x' * cov_classical * x));   
    sharpe_ratio_robust
    sharpe_ratio_classical
       
    sharpe_array_robust_classical(d_curr) = sharpe_ratio_robust / sharpe_ratio_classical;
    disp(['Mean Sharpe Ratio: ', num2str(sharpe_array_robust_classical(d_curr))]);
    
    
    % Computing Worst-Case Sharpe Ratio
    %{
    % Robust
    w = R' * x;
    cvx_begin quiet
    cvx_solver mosek
    variables pi_opt(T, 1)
    maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
    sum(kl_div(pi_opt, eta)) <= d;
    ones(1, T) * pi_opt == 1;
    pi_opt >= 0;
    cvx_end
    pi = pi_opt;
    
    sharpe_ratio_worst_robust = real(((R * pi - r_f * ones(n, 1))' * x) / sqrt(pi' * (w .* w) - pi' * (w * w') * pi));
    %}
    sharpe_ratio_worst_robust
    
    % Classical
    w = R' * x_classical;
    
    if d == 0
        pi = eta;
    else
        cvx_begin quiet
        cvx_solver mosek
        variables pi_opt(T, 1)
        maximize( pi_opt' * square(w) - quad_form(pi_opt, w * w') )
        sum(kl_div(pi_opt, eta)) <= d;
        ones(1, T) * pi_opt == 1;
        pi_opt >= 0;
        cvx_end
        pi = pi_opt;
        
    end
    cov_classical_worst = zeros(n, n);
    for k = 1 : T
        r = R(:, k);
        cov_classical_worst = cov_classical_worst + eta(k) * (r - R * eta) * (r - R * eta)';
    end
    cov = (cov_classical_worst + cov_classical_worst') / 2;
    cov_classical_worst = cov;
 
    sharpe_ratio_worst_classical = ((R * pi - r_f * ones(n, 1))' * x_classical) / sqrt(x_classical' * cov_classical_worst * x_classical);     
    sharpe_ratio_worst_classical
    sharpe_array_worst_robust_classical(d_curr) = sharpe_ratio_worst_robust / sharpe_ratio_worst_classical;
    
    disp(['Worst Case Sharpe Ratio: ', num2str(sharpe_array_worst_robust_classical(d_curr))]);

end


figure
subplot(2, 1, 1)
plot(d_array, sharpe_array_robust_classical);
xlabel('d');
ylabel('Robust/Classical Sharpe Ratio')
title('Robust/Classical Mean Sharpe Ratio');

subplot(2, 1, 2)
plot(d_array, sharpe_array_worst_robust_classical);
xlabel('d');
ylabel('Robust/Classical Sharpe Ratio')
title('Robust/Classical Worst Case Sharpe Ratio');