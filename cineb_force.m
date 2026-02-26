function [neb_f,f] = cineb_force(points,tan,spring,values,i)
% Calculating CI-NEB force
fix_index = fix_atoms;
grad = -values; % true force
grad(1,:) = zeros(1,size(points,2)); % fix reactant
grad(size(points,1),:) = zeros(1,size(points,2)); % fix product
grad(:,fix_index) = zeros(size(points,1),size(fix_index,2)); 
for index = 1:length(fix_index)
        f_ind = fix_index(index);
        grad(:,f_ind) = zeros(size(points,1),1); % initializing
end

for k=1:size(points,1)
    t=tan(k,:)';
    s=spring(k,:)';
    g = grad(k,:)'; 
    force_parallel = (dot(g,t)/((norm(t)^2)))*t;
    if k == i
        force_L = g - force_parallel;
        f(k)=norm(force_L);
        S(:,k)= g - (2*force_parallel); %-g + (2*force_parallel);
    else
        force_L = g - force_parallel;
        f(k)=norm(force_L);
        S(:,k)=(force_L)+(s); % perpendicular true force and parallel spring force
    end
    neb_f(k,:)= S(:,k)';
end

end
