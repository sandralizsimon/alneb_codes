function [alpha, optimizedParams, K, exitflag] = get_conveged_surface(Xval,yval,optimizedParams,dimension)

% Extract relevant columns
[X_ind,Y_ind] = extract_index;
values = yval(:,Y_ind); 
ndata = size(values,1);
X = Xval(:,X_ind); 
num = size(values,1)*size(values,2);
y = reshape(values,[num,1]);

pot_val = y(1:ndata);
grad_val = y(ndata+1:end);
y_values = [pot_val;grad_val];

exitflag = 10;

% Extract optimized parameters
lengthScale = optimizedParams(1)*ones(1,dimension);

row_m_sq = optimizedParams(end-2);
row_c_sq = optimizedParams(end-1);
row_e_sq = optimizedParams(end);

% Compute the kernel matrix K
K = rbf_kernel_extended_star(X, X, lengthScale, row_m_sq, row_c_sq, row_e_sq, 1);
K_inv = pinv(K);

% Compute alpha
alpha = K_inv * y_values;

end

