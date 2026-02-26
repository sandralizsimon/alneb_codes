function force_parallel=spring_parallel(points,tan,sp_const)
% Calculating parallel component of spring force
k=sp_const;
force = zeros(size(points));
force_parallel = zeros(size(points));
for i = 2 : size(points,1)-1
    force(i,:)=k*((points(i+1,:)-points(i,:)) - (points(i,:)-points(i-1,:)));
    t=tan(i,:);
    f=force(i,:);
    force_parallel(i,:) = (dot(f,t)/(norm(t)^2))*t;
end
end
