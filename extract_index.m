function [X_ind,Y_ind,len] = extract_index
% Extracting the position and values index needed for GPR learning
X_ind = 1:189; % position index
Y_ind = 1:190; % values index
len = 343*3;
end
