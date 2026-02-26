function [Ynew_inter,interpolated_points,covariance] = finer_predict(X,points,values,alpha,optimizedParams,K)
% Predict value and covariance at the finer grid 
interpolated_points = finerdist(points); % creating finer grid
[Ynew_inter,ind,covariance] = GPR_predict_pinv(X,interpolated_points,values,alpha,optimizedParams,K);
Ynew_inter = Ynew_inter*10000; %in kJ
end