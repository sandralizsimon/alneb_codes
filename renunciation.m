%% Renunciation

%% Parameters for EMT potential 
De = 0.7102;    % (eV)
a  = 1.6047;    % (1/A)
re = 2.897;     % (A)
rc = 9.5;       % (A)

%% Read reactant structure from POSCAR
poscar_data = readPOSCAR('POSCAR_initial.1');

box          = poscar_data.lattice_vectors;   % Simulation cell vectors
len          = poscar_data.atom_counts;       % Number of atoms
fixatoms_old = poscar_data.flags;             % Fixed/free atom flags
positions1  = poscar_data.positions;          % Atomic coordinates

%% Compute energy and forces at reactant
[energy1, forces1] = find_forces_emt(positions1, box, De, a, re, rc, fixatoms_old, 0);

%% Generate updated POSCAR positions based on forces
[new_positions1, fixatoms] = create_new_poscar(len, forces1, fixatoms_old, positions1);

%% Load data from Exploration-Exploitation loop
Xval = load("finaldata.txt");      % Atomic configurations 
yval = load("finalvalues.txt");    % Energies + forces

values = yval;
values(:,2:end) = -1*yval(:,2:end); 
% Convert forces → gradients (−F = ∇E)

%% Initial Transition State (TS) guess from from Exploration-Exploitation loop
TS_guess = Xval(end-2,:); 

%% High-energy (HE) images around TS guess
HE_images = [Xval(end-1,:); Xval(end-2,:); Xval(end,:)];

%% Compute distances of all points from TS guess
distances = pdist2(Xval, TS_guess);

%% Create hypersphere around TS
hemi_points = [];
hemi_values = [];

dist_limit = max(distances(end-2:end)); % Radius defined by HE images

for ii = 1:length(distances)
    if distances(ii) <= dist_limit
        hemi_points = [hemi_points; Xval(ii,:)];
        hemi_values = [hemi_values; values(ii,:)];
    end
end
disp(['Number of data points in hemisphere: ', num2str(size(hemi_points,1))])

%% GPR hyperparameters
dimension = 189;        % Input dimension (3N for N atoms)
initialLengthScale = 0.8; % Lengthscale
row_m_sq = 1;           % variance
row_c_sq = 0.005;       % Energy noise 
row_e_sq = 0.0005;      % Force noise 

initialParams = [initialLengthScale, row_m_sq, row_c_sq, row_e_sq];

%% Transition-state convergence settings
max_TS_gradient = 10;   % Initial large gradient
nneval = 200;           % Energy/force evaluation counter

%% Extract gradient at initial TS estimate
TS_gradient = hemi_values(end-2,2:end) ./ 96.4853; 
% Convert forces from kJ/mol/A to eV/Å

prev_max_grad = abs(max(TS_gradient, [], "all", "ComparisonMethod", "abs"));

%% Initialize TS data structure
TS = 0;                 % TS convergence flag

ts_data.points        = HE_images(2,:);   % TS coordinates
ts_data.gradient_max  = prev_max_grad;    % Max gradient norm
ts_data.values        = hemi_values(end-2,:);
cipoints              = HE_images;
ts_data.mep           = cipoints;          % Initial MEP
ts_data.potval        = [];                % Predicted energies

nr_iter   = 1;
iteration = 1;
optimizedParams = initialParams;

%% Display initial renunciation summary
disp('=== Renunciation Phase Started ===')
disp(['Initial TS max gradient: ', num2str(prev_max_grad), ' eV/A (target < 0.02)'])
disp(['Hemisphere radius: ', num2str(dist_limit), ' A'])
disp(['Initial hemisphere data points: ', num2str(size(hemi_points,1))])
disp(['Initial TS energy: ', num2str(hemi_values(end-2,1)), ' kJ/mol'])

%% Main TS optimization loop
while TS == 0
    disp(['--- Renunciation Iteration ', num2str(iteration), ' ---'])
    tic
    disp('Learning GPR model...')
    % Train GPR model using hemispherical data
    [alpha, optiparam, K, exitflag] = GPR_alpha_pinv_bayes(hemi_points, hemi_values, optimizedParams, dimension);
    toc

    % Handle failed hyperparameter optimization
    if exitflag == -1
        [alpha, optimizedParams, K, exitflag] = GPR_alpha_pinv_bayes(hemi_points, hemi_values, initialParams, dimension);
    else
        optimizedParams = optiparam;
    end

    %% Predict energies on climbing-image points
    cipoints = HE_images;
    ciimg = 2;                    % Index of climbing image

    %% Perform CI-NEB optimization
    [potval, mep, iter, neb_flag, max_dist] = GPR_ci_neb_opt_pinv( hemi_points, hemi_values, cipoints,alpha, optimizedParams,K, 4000, 0.2, ciimg);

    % Fallback with smaller step size if NEB fails
    if iter == 1000 || neb_flag == 0
        [potval, mep, iter, neb_flag, max_dist] = GPR_ci_neb_opt_pinv(hemi_points, hemi_values, cipoints, alpha, optimizedParams, K, 50, 0.1, ciimg);
    end

    %% Extract updated TS geometry
    x_ts = mep(ciimg,:);
    iteration = iteration + 1;

    %% Evaluate true energy and forces at TS
    bond_pos = points_to_pos(x_ts);
    nneval = nneval + 1;
    [energyn, forcesn] = find_forces_emt(bond_pos, box, De, a, re, rc, fixatoms, nneval);

    current_value_ts = [energyn, -pos_to_points(forcesn)];

    %% Compute TS gradient norm
    TS_gradient = current_value_ts(2:end) ./ 96.4853;
    max_TS_gradient = abs(max(TS_gradient, [], "all", "ComparisonMethod", "abs"));

    disp(['  Max TS gradient: ', num2str(max_TS_gradient), ' eV/A (target < 0.02)'])
    disp(['  TS energy: ', num2str(current_value_ts(1)), ' kJ/mol'])

    %% Store TS history
    ts_data.gradient_max = [ts_data.gradient_max; max_TS_gradient];
    ts_data.points       = [ts_data.points; x_ts];
    ts_data.values       = [ts_data.values; current_value_ts];
    ts_data.mep          = [ts_data.mep; mep];
    ts_data.potval       = [ts_data.potval; potval];

    save('ts_data_final.mat','ts_data')

    %% Check convergence
    if max_TS_gradient < 0.02
        disp('=== TS Converged! ===')
        disp(['  Final max gradient: ', num2str(max_TS_gradient), ' eV/A'])
        disp(['  Final TS energy: ', num2str(current_value_ts(1)), ' kJ/mol'])
	    disp([' Training points added in renunciation phase: ', num2str(iteration-1)])
        TS = 1;   % TS converged
        break;
    else
        % Add new TS point to training set
        hemi_points = [hemi_points; x_ts];
        hemi_values = [hemi_values; current_value_ts];
        disp('Data point added to hemisphere')
        % Update HE image if gradient improves
        if max_TS_gradient < prev_max_grad
            HE_images(2,:) = x_ts;
            prev_max_grad = max_TS_gradient;
            disp('  HE_images updated (gradient improved)')
        end
    end
end
