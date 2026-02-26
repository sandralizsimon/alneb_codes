function [neb_f,f] = neb_force(points,tan,spring,values)
% Calculating NEB force
fix_index = fix_atoms;
grad = -values; % true force
grad(1,:) = zeros(1,size(points,2)); % fixing reactant
grad(size(points,1),:) = zeros(1,size(points,2)); % fixing product
grad(:,fix_index) = zeros(size(points,1),size(fix_index,2)); % initializing
for k=1:size(points,1)
    t=tan(k,:)';
    s=spring(k,:)';
    g = grad(k,:)';
    force_parallel = (dot(g,t)/(norm(t)^2))*t;
    force_L = g - force_parallel;
    f(k)=norm(force_L);
    S(:,k)=(force_L)+(s); % perpendicular true force and parallel spring force
    neb_f(k,:)= S(:,k)';
end
end