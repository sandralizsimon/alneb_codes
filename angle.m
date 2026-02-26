function angle_in_degrees=angle(v1,v2)

% Calculate the dot product of the two arrays
dot_product = dot(v1(:), v2(:));

% Calculate the norm (magnitude) of each array
norm_A = norm(v1(:));
norm_B = norm(v2(:));

% Calculate the angle in radians
angle_in_radians = acos(dot_product / (norm_A * norm_B));

% Convert the angle to degrees
angle_in_degrees = rad2deg(angle_in_radians);

end
