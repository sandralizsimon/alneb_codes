function [alpha, optimizedParams, K, exitflag] = GPR_alpha_pinv_bayes(Xval,yval,initialParams,dimension)
% Optimising hyperparameters
% Extract relevant columns
[X_ind,Y_ind,len] = extract_index;
values = yval(:,Y_ind); % For 2D
X = Xval(:,X_ind); % For 2D
num = size(values,1)*size(values,2);
y = reshape(values,[num,1]);

objFcn = @(params) negLogLikelihood_pinv_bayesian(params,X,y,dimension);

% Define constraints to ensure positivity of length scale and noise variance
nonlcon = @(params) nonlconFunction(params);

try
    lb = [0.2,0.2, 0.005, 0.0005]; % Lower bounds for each parameter
    ub = [1.5, 1, 0.05 ,0.005]; % Upper bounds for each parameter
    % Create parallel pool if it's not already created
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Set options with parallel computing enabled
    options = optimoptions('fmincon', 'UseParallel', true, 'Algorithm','interior-point','Display', 'iter');
    [optimizedParams,~,exitflag, ~] = fmincon(objFcn, initialParams,[],[],[],[], lb, ub,[], options);
    
    % Check if any of the optimized parameters are Inf, NaN, or imaginary
    if any(isinf(optimizedParams)) || any(isnan(optimizedParams)) || ~isreal(optimizedParams)
        error('Optimization resulted in Inf, NaN, or imaginary values.');
    end

catch ME
    % Handle any errors that occur during optimization
    fprintf('An error occurred during optimization: %s\n', ME.message);
    optimizedParams = zeros(size(initialParams));
    exitflag = -1;
end

% Extract optimized parameters
lengthScale = optimizedParams(1)*ones(1,dimension);
row_m_sq = optimizedParams(end-2);
row_c_sq = optimizedParams(end-1);
row_e_sq = optimizedParams(end);

% Compute the kernel matrix K
n = length(X);
K = rbf_kernel_extended_star(X, X, lengthScale, row_m_sq, row_c_sq, row_e_sq, 1);
K_inv = pinv(K);

% Compute alpha
alpha = K_inv * y;

end
