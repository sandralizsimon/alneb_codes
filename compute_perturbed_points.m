function [perturbed_point1,perturbed_point2] = compute_perturbed_points(point1,point2)
dist = 0.2;

% Perturb point1 along the vector connecting point1 and perturbed_midpoint
perturbed_point1 = point1 + dist * (point2 - point1) / norm(point2 - point1);

% Perturb point2 along the vector connecting point2 and perturbed_midpoint
perturbed_point2 = point2 + dist * (point1 - point2) / norm(point1- point2);

end