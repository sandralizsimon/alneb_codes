function new_value = find_high_energy(image_values,neb_final,X,values,alpha,optimizedParams,K)
% Finding the maximum predicted energy configuration
pot_energy = [-inf;image_values(2:size(neb_final,1)-1,1);-inf]; % only checking 2nd to 2nd last image
[max_potential,index] = max(pot_energy); 
max1 = neb_final(index, :);
max0 = neb_final(index-1, :);
max2 = neb_final(index+1, :);

% number of grid points in each segment to check for maximum
nx = 10;

for i = 1:nx
    t = (i - 1) / (nx - 1);  % Parameter 't' ranges from 0 to 1
    interpolatedPoints1(i, :) = (1 - t) * max0 + t * max1;
    interpolatedPoints2(i, :) = (1 - t) * max1 + t * max2;
end
images = [interpolatedPoints1;interpolatedPoints2(2:end,:)];

[Yinterpolated,ind,cov_images] = GPR_predict_pinv(X,images,values,alpha,optimizedParams,K);
image_v = Yinterpolated(:,1);
pot_energy = [-inf;image_v(2:size(images,1)-1,1);-inf]; % only checking 2nd to 2nd last image
[max_val,max_index] = max(pot_energy(:,1));
new_value = images(max_index,:);

end