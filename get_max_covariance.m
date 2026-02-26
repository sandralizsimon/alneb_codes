function [new_value,max_idx,covariance] = get_max_covariance(points,covariance,data)
% Find the configuration with maximum  predicted covariance
flag = 0;
while flag == 0
    [mval,ival] = max(covariance);
    hval = find(covariance==mval);
    randomIndex = randi(numel(hval));
    randomCovInd = hval(randomIndex);
    new_value = points(randomCovInd,:);
    max_idx = randomCovInd;

    if any(pdist2(new_value,data) <= 0.001)
        flag = 0;
        covariance(max_idx) = -inf;
        if all(covariance == -inf)
            new_value = [];
            flag = 3;
        end
    else
        flag = 1;
    end
end
end