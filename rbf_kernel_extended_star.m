% Define the RBF kernel function
function final_kernal = rbf_kernel_extended_star(x1, x2,ld,row_m_sq,row_c_sq,row_e_sq,ind)
lengthscale = ld.^2;
a = row_m_sq;
b = row_c_sq.^2;
c = row_e_sq.^2;

num_samples_x1 = size(x1, 1);
num_samples_x2 = size(x2, 1);
K_x_x = zeros(num_samples_x1, num_samples_x2);

for i = 1:num_samples_x1
    for j = 1:num_samples_x2
        sum_squared_diff = 0;
        for k = 1:size(x1, 2)
            sum_squared_diff = sum_squared_diff + (((x1(i, k) - x2(j, k))^2)/lengthscale(k));
        end
        K_x_x(i, j) = exp(-0.5 * sum_squared_diff);
    end
end

% cov (f(x1), f'(x1))
num_samples_x1 = size(x1, 1);
num_samples_x2 = size(x2, 1);
num_dimensions = size(x1, 2);

K_x_dx = zeros(num_samples_x1, num_samples_x2 * num_dimensions);

for d = 1:num_dimensions
    result_matrix2 = zeros(num_samples_x1, num_samples_x2);

    for i = 1:num_samples_x1
        for j = 1:num_samples_x2
            result_matrix2(i, j) = a*K_x_x(i,j)* (1 / lengthscale(d)) *(x1(i, d) - x2(j, d));
        end
    end
    K_x_dx(:, (d - 1) * num_samples_x2 + 1:d * num_samples_x2) = result_matrix2;
end

% cov (f'(x1), f(x2))
num_samples_x1 = size(x1, 1);
num_samples_x2 = size(x2, 1);
num_dimensions = size(x1, 2);

K_dx_x = zeros(num_samples_x1 * num_dimensions, num_samples_x2);

for d = 1:num_dimensions
    result_matrix3 = zeros(num_samples_x1, num_samples_x2);

    for i = 1:num_samples_x1
        for j = 1:num_samples_x2
            result_matrix3(i, j) = a*K_x_x(i,j)* (1 / lengthscale(d)) * (-(x1(i, d) - x2(j, d)));
        end
    end
    K_dx_x((d - 1) * num_samples_x1 + 1 : d * num_samples_x1, :) = result_matrix3;
end

% cov (f'(x1), f'(x2))

K_dx_dx = [];
Kxx = [];
for d1 = 1:num_dimensions
    for d2 = 1:num_dimensions
        for i = 1:num_samples_x1
            for j = 1:num_samples_x2
                sum_squared_diff = 0;
                if d1 == d2
                    val = 1;
                else
                    val = 0;
                end
                result_matrix4(i, j) = a*K_x_x(i,j)*(1/lengthscale(d1))*(val - ((x1(i, d1) - x2(j, d1))*(x1(i, d2) - x2(j, d2)))/lengthscale(d2));
            end
        end
        Kxx = [Kxx , result_matrix4];
    end
    K_dx_dx = [K_dx_dx;Kxx];
    Kxx = [];
end

K_x_x = (a*K_x_x);
K_dx_dx = (K_dx_dx);
if ind == 1
    final_kernal = [[(K_x_x+ (b)*eye(size(K_x_x))),K_x_dx];[K_dx_x,(K_dx_dx + (c)*eye(size(K_dx_dx)))]];
elseif ind == 2
    final_kernal = [K_x_x,K_x_dx];
elseif ind == 3
    final_kernal = [K_dx_x,K_dx_dx];
elseif ind == 4
    final_kernal = K_x_x;
elseif ind == 5
    final_kernal = K_dx_x;
end

end
