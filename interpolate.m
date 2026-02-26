function points = interpolate(point1, point2)

% Linear interpolation
n = 7;
interpolatedPoints1 = zeros(n, length(point1));
m1 = point1;
m2 = point2;

for i = 1:n
    t = (i - 1) / (n - 1);  % Parameter 't' ranges from 0 to 1
    interpolatedPoints1(i, :) = (1 - t) * m1 + t * m2;
end

points = interpolatedPoints1; % points in Angstroms  

end
