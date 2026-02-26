function interpolated_points = finerdist(points)
% Number of points
num_points = size(points, 1);

% Number of dimensions
num_dimensions = size(points, 2);

% Number of interpolated points between each pair
num_interpolated_points = 2;

% Initialize interpolated points matrix
interpolated_points = zeros((num_points-1)*(num_interpolated_points+1) + 1, num_dimensions);

% Interpolate between each pair of neighboring points
for i = 1:num_points-1
    % Extract neighboring points
    p1 = points(i,:);
    p2 = points(i+1,:);

    % Interpolate between p1 and p2
    interpolated = interp1([0, 1], [p1; p2], linspace(0, 1, num_interpolated_points+2));

    % Assign interpolated points to the result matrix
    start_idx = (i-1)*(num_interpolated_points+1) + 1;
    end_idx = start_idx + num_interpolated_points + 1;
    interpolated_points(start_idx:end_idx, :) = interpolated;
end

% The last point is the same as the last point in the original data
interpolated_points(end,:) = points(end,:);

end
