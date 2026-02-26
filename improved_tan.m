function tan = improved_tan(points,P)
% Improved tangent NEB
c=zeros(size(points));
R=points;
n=size(points,1);

c(1,:)=(R(2,:)-R(1,:))/norm(R(2,:)-R(1,:));
c(n,:)=(R(n,:)-R(n-1,:))/norm(R(n,:)-R(n-1,:));

for i=2:n-1
    if P(i+1)>P(i) && P(i)>P(i-1)
        t = R(i+1,:) - R(i,:);
        c(i,:)=t/norm(t); 
    elseif P(i+1)<P(i) && P(i)<P(i-1)
        t = R(i,:) - R(i-1,:);
        c(i,:)=t/norm(t);
    elseif (P(i+1)<P(i) && P(i)>P(i-1)) || (P(i+1)>P(i) && P(i)<P(i-1))
        
        v_max = max(norm(P(i+1) - P(i)),norm(P(i-1) - P(i)));
        v_min = min(norm(P(i+1) - P(i)),norm(P(i-1) - P(i)));
        b1 = R(i+1,:) - R(i,:);
        b2 = R(i,:) - R(i-1,:);
        if P(i+1)>P(i-1)
            t=(b1*v_max)+(b2*v_min);
            c(i,:)=t/norm(t);
        elseif P(i+1)<P(i-1)
            t=(b1*v_min)+(b2*v_max);
            c(i,:)=t/norm(t);
        else
            t=(b1*v_min)+(b2*v_max);
            c(i,:)=t/norm(t);
        end                  
        
    end
end

tan = c;

end