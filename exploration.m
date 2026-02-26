%% Exploration
tic
%% Parameters for EMT potential (values for Pt)
De = 0.7102;          % Depth of potential well (eV)
a = 1.6047;           % Width parameter (1/A)
re = 2.897;           % Equilibrium bond length (A)
rc = 9.5;            % Cutoff radius (A)

disp('Reading initial and final geometries');
%% Read Initial Geometry
poscar_data = readPOSCAR('POSCAR_initial.1');
box = poscar_data.lattice_vectors;
len = poscar_data.atom_counts;
fixatoms_old = poscar_data.flags;
positions1 = poscar_data.positions;
% Calculate forces and relax the structure using EMT
[energy1, forces1] = find_forces_emt(positions1, box, De, a, re, rc, fixatoms_old, 0);
[new_positions1, fixatoms] = create_new_poscar(len, forces1, fixatoms_old, positions1);
[energy1, forces1] = find_forces_emt(new_positions1, box, De, a, re, rc, fixatoms, 1);
point1 = pos_to_points(new_positions1); % Flatten coordinates to a vector

%% Read Final Geometry
poscar_data = readPOSCAR('POSCAR_final.1');
positions2 = poscar_data.positions;
[energy2, forces2] = find_forces_emt(positions2, box, De, a, re, rc, fixatoms_old, 100);
[new_positions2, fixatoms] = create_new_poscar(len, forces2, fixatoms_old, positions2);
[energy2, forces2] = find_forces_emt(new_positions2, box, De, a, re, rc, fixatoms, 2);
point2 = pos_to_points(new_positions2);

% Format initial training data (Energy + Force vector)
pos_trial = pos_to_points(forces1);
gaussian_values = [energy1, zeros(size(pos_trial)); energy2, zeros(size(pos_trial))];

%% Compute Perturbed minima 
disp("Computing perturbed minima");
[perturbed_point1, perturbed_point2] = compute_perturbed_points(point1, point2);
perturbed_pos1 = points_to_pos(perturbed_point1);
perturbed_pos2 = points_to_pos(perturbed_point2);

%% Starting Dataset (X = Positions, gaussian_values = E & F)
X = [point1; point2; perturbed_point1; perturbed_point2];
[energy3, forces3] = find_forces_emt(perturbed_pos1, box, De, a, re, rc, fixatoms, 3);
[energy4, forces4] = find_forces_emt(perturbed_pos2, box, De, a, re, rc, fixatoms, 4);
gaussian_values(3,:) = [energy3, pos_to_points(forces3)];
gaussian_values(4,:) = [energy4, pos_to_points(forces4)];
eval = 4; % Total number of true force evaluations performed

% Save initial datasets
writematrix(gaussian_values, 'finalvalues.txt', 'Delimiter', 'tab');
writematrix(X, 'finaldata.txt', 'Delimiter', 'tab');

%% Hyperparameter Initialization for GPR
dimension = 189;             % DOF
initialLengthScale = 0.8;    % RBF kernel lengthscale
row_m_sq = 1;                % variance
row_c_sq = 0.001;            % Energy noise 
row_e_sq = 0.001;            % Force noise 
initialParams = [initialLengthScale, row_m_sq, row_c_sq, row_e_sq];

%% Variable Initialization for the Main Loop
start = 1; 
learn = 1;                   % Toggle for training
iteration = 0; 
param = [];                  % Storage for optimized hyperparameters
exits = [];                  % Storage for optimization exit flags
CINEB = 0;                   % Flag for Climbing Image NEB phase
neb_flag = 0;                % Flag for NEB convergence
d_max = 10;                  % Convergence metric (MEP movement) (initializing)
d_min = 10;                  % Min distance to existing training data (initializing)
max_iterations = 500;        % Max steps for the NEB optimizer
points = interpolate(point1, point2); % Linear interpolation between minima
initial_mep = points; 
neb_previous = points;
neb_final = points;
rmax_limit = 2.5;            % Max allowed deviation from initial path
dmin_limit = 0.1;               % Min distance to existing training data
dmax_limit = 0.2;              % Convergence threshold for MEP movement
initial_param_values = initialParams;

% Counter variables for error handling and stalling
rcount = 0; ncount = 0; del_count = 0; count = 0;
counter_limit = 3; counter_sum = 5;
new_neb_flag = 0;
skip_learning = 0;
count_perturb = 0;
perturb_iter = 50;
status = 0;
neb_paths = [];

%% Output File Header Initialization
array_out = {'iteration','exitflag','neb_flag','perturbation','neb_iterations','rmax','dmax', ...
    'dmin','mep','predicted_pot'};

rng("shuffle"); % random seed 

%% MAIN  LOOP
while d_max > dmax_limit || iteration <= 5  ||  CINEB == 0
    s = ['iteration','_',num2str(iteration)];
    disp(s);
    
    % Evaluate true energy and forces
    if (iteration > 0) && learn == 1
        bond_pos = points_to_pos(X(end,:));
        eval = eval + 1;
        [energyn, forcesn] = find_forces_emt(bond_pos, box, De, a, re, rc, fixatoms, eval);
        current_value = [energyn, pos_to_points(forcesn)];
        gaussian_values = [gaussian_values; current_value];
        writematrix(gaussian_values, 'finalvalues.txt', 'Delimiter', 'tab');
    end
    learn = 1;

    if status <= 0 
        % GPR Hyperparameter Optimization
        figure('visible', 'off');
        disp("Start optimization")
        values = gaussian_values;
        values(:,2:end) = -1*gaussian_values(:,2:end); % Convert forces to gradients for GPR
        
        if skip_learning == 1
            [alpha_value, finalParams, K_matrix, exitflag] = get_conveged_surface(X, values, optimizedParams, dimension);
            skip_learning = 0;
        else
            [alpha_value, finalParams, K_matrix, exitflag] = GPR_alpha_pinv_bayes(X, values, initial_param_values, dimension);
        end
        disp("Optimization ended")
        param = [param; finalParams];
        exits = [exits; exitflag];

        if exitflag >= 0 % Successful GPR Optimization
            disp("Optimised Parameters Found")
            disp(finalParams)
            optimizedParams = finalParams;
            alpha = alpha_value;
            K = K_matrix;
            
            % NEB optimization on the GPR surrogate surface
            if new_neb_flag == 0 % Standard NEB on surrogate surface
                perturb = 0;
                [potval, mep, iter, neb_flag] = GPR_neb_opt_pinv(X, values, initial_mep, alpha, optimizedParams, K, max_iterations);
            elseif new_neb_flag == 1
                disp("perturbed mep") % NEB with perturbation
                perturb = 1;
                [potval, mep, iter, neb_flag] = GPR_neb_opt_pinv(X, values, initial_mep, alpha, optimizedParams, K, perturb_iter);
                new_neb_flag = 0;
            end

            % Step 4: Analyze NEB Results
            if neb_flag == 1 % NEB on surrogate surface converged
                disp('NEB converged')
                rmax = distance(mep, points); % Distance between MEP and initial guess
                
                if (rmax <= rmax_limit) 
                    disp("rmax within limit")
                    neb_final = mep;
                    neb_paths = [neb_paths; [iteration, (potval(:,1))']];
                    
                    % Find configuration with highest uncertainty (covariance) to sample next
                    [potval_inter, finer_mep_final, covariance_final] = finer_predict(X, neb_final, values, alpha, optimizedParams, K); % Finer-grid
                    [new_value, max_idx, covariance_final] = get_max_covariance(finer_mep_final, covariance_final, X);
                    
                    % Reset counters 
                    count = 0; rcount = 0; ncount = 0; del_count = 0; count_perturb = 0;
                    initial_mep = neb_final;
                    d_max = distance(neb_previous, neb_final); % Track MEP movement between iterations
                    neb_previous = neb_final;
                    initial_param_values = optimizedParams;
                else
                    disp('rmax beyond limit') 
                    [potval_inter, finer_mep_final, covariance_final] = finer_predict(X, neb_final, values, alpha, optimizedParams, K);
                    [new_value, max_idx, covariance_final] = get_max_covariance(finer_mep_final, covariance_final, X);
                    rcount = rcount + 1;
                    if rcount >= counter_limit
                        disp('Re-initialize')
                        initial_mep = points;
                        initial_param_values = initialParams;
                        if rcount + ncount >= counter_sum
                            new_neb_flag = 1;
                        end
                    end
                end
            else
                % NEB failed to converge on the surrogate
                disp('NEB not converged')
                [potval_inter, finer_mep_final, covariance_final] = finer_predict(X, neb_final, values, alpha, optimizedParams, K);
                [new_value, max_idx, covariance_final] = get_max_covariance(finer_mep_final, covariance_final, X);
                ncount = ncount + 1;
                if ncount >= counter_limit
                    disp('Re-initialize')
                    initial_mep = points;
                    initial_param_values = initialParams;
                    if ncount + rcount >= counter_sum
                        new_neb_flag = 1;
                    end
                end
            end
        else
            % Optimization of hyperparameters failed
            disp('Optimised Parameters not found')
            [new_value, max_idx, covariance_final] = get_max_covariance(finer_mep_final,covariance_final, X);
            count = count + 1;
            if count >= counter_limit
                disp('Re-initialize')
                initial_mep = points;
                initial_param_values = initialParams;
                if count >= counter_limit
                    new_neb_flag = 1; 
                    skip_learning = 1;
                end
            end

        end
    else
        % Calculation error handling
        X(end,:) = [];
        [new_value, max_idx, covariance_final] = get_max_covariance(finer_mep_final,covariance_final, X);
        del_count = del_count + 1;
        if del_count >= counter_limit
            disp('Re-initialize')
            initial_mep = points;
            if del_count >= counter_sum
                new_neb_flag = 1;
                count_perturb = count_perturb + 1;
            end
        end
    end

    % Logging and Plotting
    newRow1 = {iteration, exitflag, neb_flag, perturb, iter, rmax, d_max, d_min, mep, potval};
    array_out(end+1,:) = newRow1;
    save('array_out.mat', 'array_out');
    X = [X; new_value]; % Add new configuration to the training set
    iteration = iteration + 1;
    writematrix(X, 'finaldata.txt', 'Delimiter', 'tab');
    
    %% Checking convergence
    if d_max < dmax_limit && iteration > 5 && exitflag >= 0
        disp("2nd Phase - Exploitation")
        X(end,:) = []; % remove the last selection
        disp("Highest energy point")
        new_value = find_high_energy(potval(:,1), neb_final, X, values, alpha, optimizedParams, K);
        d_min = min(pdist2(new_value, X));
        
        % Starting Exploitation
        [HE_images, d_max, X, values, gaussian_values, alpha, optimizedParams, K, iteration, eval,neb_paths, array_out] = exploitation(new_value, X, values, gaussian_values, alpha, optimizedParams, K, eval, d_min, d_max, initialParams, initial_param_values, max_iterations, points, rmax_limit, dmin_limit, counter_limit, counter_sum, dmax_limit, neb_final, iteration, max_idx, dimension, neb_paths, potval, box, De, a, re, rc, fixatoms, array_out);
      
        status_list = [];
        HE_energies = inf*ones(size(HE_images,1), size(HE_images,2)+1);
        maxforce_at_ts = 0;
       
        % Bracketing Check: Ensure the TS candidate is actually a maximum
        display("Bracketing Check");
        for k = [2, 1, 3]
            if maxforce_at_ts < 0.25
                new_value = HE_images(k,:);
                display(['Current HE_images(', num2str(k), ',:) = ', num2str(new_value(1:5))]);
                X = [X; new_value];
                writematrix(X, 'finaldata.txt', 'Delimiter', 'tab');
                eval = eval + 1;
                bond_pos = points_to_pos(X(end,:));
                [energyn, forcesn] = find_forces_emt(bond_pos, box, De, a, re, rc, fixatoms, eval);
                current_value = [energyn, pos_to_points(forcesn)];
                gaussian_values = [gaussian_values; current_value];
                writematrix(gaussian_values, 'finalvalues.txt', 'Delimiter', 'tab');
                HE_energies(k,:) = current_value;
                display(['Current HE_energies(', num2str(k), ',1) = ', num2str(HE_energies(k,1))]);
                status_list = [status_list; status];
                force_at_ts = HE_energies(2, 2:end)./96.4853;
                maxforce_at_ts = abs(max(force_at_ts, [], "all", 'ComparisonMethod', "abs"));
                display(['Current maxforce_at_ts = ', num2str(maxforce_at_ts)]);
            end
        end

        % Final Convergence Check: If TS image is higher than its neighbors, we switch to CI-NEB
        if (HE_energies(1,1) < HE_energies(2,1)) && (HE_energies(3,1) < HE_energies(2,1)) && all(status_list) == 0
            CINEB = 1;
            disp(['Number of training points excluding minima: ', num2str(size(gaussian_values,1)-2)]);
            disp("3rd Phase - Renunciation")
            renunciation; % Start renunciation
            break;
        else
            CINEB = 0;
            initial_param_values = optimizedParams;
            learn = 0;
            d_min = 10;
            d_max(end) = 10;
            skip_learning = 0;
        end
    end
    
end