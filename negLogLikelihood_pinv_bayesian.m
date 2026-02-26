function [negLL,K,K_inv,alpha] = negLogLikelihood_pinv_bayesian(params, X, y,dimension)
ld = params(1)*ones(1,dimension);
row_m_sq = params(end-2);
row_c_sq = params(end-1);
row_e_sq = params(end);


detK  = 0;
i =0;
while  detK == 0 
    if i ~= 0
        row_c_sq = row_c_sq*10;
        row_e_sq = row_e_sq*10;
    end

    K = rbf_kernel_extended_star(X, X,ld,row_m_sq,row_c_sq,row_e_sq,1); %+ lambda * eye(n);
    detK = det(K);
    i = i+1;
end
detK = det(K);
log_det_K = log(detK);

if det(K) == inf
    L = chol(K, 'lower');  % Cholesky decomposition: K = L * L'
    log_det_K = 2 * sum(log(diag(L)));  % log(det(K)) = 2 * sum(log(diag(L)))
end

n = size(K,1);

% Calculate Kernel inverse
K_inv = pinv(K);

% Compute alpha
alpha = K_inv * y;

save('K_not_conv.mat','K');

% Negative log-likelihood
negLL = 0.5 * ((y' * alpha) + log_det_K + (n * log(2 * pi)));


end
