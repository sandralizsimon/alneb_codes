function [new_val,ind,covariance,K_star] = GPR_predict_pinv(X,Xnew,yval,alpha,optimizedParams,K)
% Calculate GPR predicted energies and forces for running NEB on model surface
[X_ind,Y_ind,len] = extract_index;
lengthScale =optimizedParams(1)*ones(1,size(X_ind,2));
row_m_sq = optimizedParams(end-2);
row_c_sq = optimizedParams(end-1);
row_e_sq = optimizedParams(end);
X = X(:,X_ind);
ndata = size(X,1);
Xnew = Xnew(:,X_ind);
y = yval(:,Y_ind);

K_star = rbf_kernel_extended_star(Xnew, X, lengthScale, row_m_sq, row_c_sq, row_e_sq, 2);
fnew_normalized = K_star * alpha;
fnew = fnew_normalized;

K_starstar = rbf_kernel_extended_star(Xnew,Xnew,lengthScale,row_m_sq,row_c_sq,row_e_sq,4);
Cov = K_starstar - (K_star * pinv(K) * K_star');
covariance = diag(Cov);
covariance(1) = -inf;
covariance(end) = -inf;
[maximum_var,ind] = max(covariance);

dK_star = rbf_kernel_extended_star(Xnew,X, lengthScale,row_m_sq,row_c_sq,row_e_sq,3);
dfnew = dK_star * alpha;
dfnew = reshape(dfnew,size(Xnew));
dfnew_final = zeros(size(Xnew,1),len);
dfnew_final(1:size(Xnew,1),X_ind) = dfnew; 

new_val = 0.0001*[fnew,dfnew_final];  % kJ/mol to (g A^2 /fs^2)/mol

end
