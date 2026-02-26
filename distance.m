function maxd = distance(oldneb,newneb)
% Calulating maximum atomic displacement between 2 neb pathways
[images,dim] = size(newneb);
d = [];
summ = zeros(1,3);
for i= 1:images
    k = 1;
    for j = 1:3:dim
        d1 = oldneb(i,j:j+2) - newneb(i,j:j+2);
        summ(k) = sqrt(sum(d1.^2));
        k = k+1;
    end
    d = [d,summ];
end
maxd = max(d);
end

