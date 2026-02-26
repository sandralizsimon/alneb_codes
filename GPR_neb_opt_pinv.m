function [Ynew,points,i,neb_flag] =  GPR_neb_opt_pinv(X,values,points,alpha,optimizedParams,K,max_iterations)
tic
disp('Starting NEB optimisation')

con_fact = 103.6484; % from gA/ft^2 mol  to eV
%% Predicting values at points using fitted GPR
[Ynew,ind,covariance] = GPR_predict_pinv(X,points,values,alpha,optimizedParams,K);
energy = Ynew(:,1);
force = Ynew(:,2:end);
neb_it = 1;

%% Tangent evaluation
tan=improved_tan(points,energy);

%% Calculating spring forces
sp_const = 5 * 0.009648; % from eV to gA/ft^2 mol  
spring = spring_parallel(points,tan,sp_const);

%% Calculating NEB forces at points
[neb_f,f] = neb_force(points,tan,spring,force);

maxforce(neb_it) = max(max(neb_f)).*con_fact; % maximum NEB force in eV/A
% Open the file for writing
fileID = fopen('maxforce_gpr.txt', 'w');
% Write initial value to file
fprintf(fileID, '%d\n', maxforce(neb_it));

% Open the file for writing
fileID1 = fopen('normforce_gpr.txt', 'w');


%% NEB evaluation using AARE

%%%AARE parameters
AARE_dt=0.1; %fs
AARE_fdec=0.5;
AARE_finc = 1.1;
AARE_dtmax = 1;
AARE_mass = ones(size(neb_f)); 

v=zeros(size(points,1),size(points,2));
grad1 = neb_f;
i=1;
k=0;
an=0;
normforce(neb_it) = norm(neb_f.*con_fact);
disp(['Initial NEB force = ', num2str(normforce(neb_it))])
% Write initial value to file
fprintf(fileID1, '%d\n', normforce(neb_it));


while  ((normforce(neb_it)> 0.05) && (i < max_iterations))
iter =i;
    if i>1
        % Fletcher Reeves
        gg=grad1;
        value=dot(grad1,gg);
        value1=dot(grad2,grad2);
        be_FR=value/value1;

        % HS
        value=dot(grad1,grad1-grad2);
        value1=dot(CG_f,grad1-grad2);
        be_HS =value/value1;

        % Displaying angles on the plot
        an=(angle(grad1,CG_f));
        ang=num2str(round(an));

        % Algorithm selection based on angles
        if an < 90
            be=be_FR;
        elseif (an > 90) && (an < 120)
            be=be_HS;
        elseif an > 120
            be=0;
        end

        % Direction Update
        CG_f= (grad1) + (be*CG_f);
    else
        CG_f = neb_f ;
    end

    v = CG_f/norm(CG_f); % velocity modification

    % for acceleration based on angles
    if an<90
        if i>1 %%% This is done to prevent accelerating in the first step...
            AARE_dt = min( AARE_dt * AARE_finc, AARE_dtmax );
        end

    elseif (an>90) && (an<120)

        AARE_dt = AARE_dt * AARE_fdec;

    end

    %%Euler integrator x_i+1 = x_i + v_i dt + 1/2 a_i dt^2
    p_new = points + (AARE_dt * v) + ((AARE_dt^2)*(neb_f)./(2*AARE_mass)); %%Euler integrator x_i+1 = x_i + v_i dt + 1/2 a_i dt^2

    %%% v_i+1 = v_i(mod) + Force/m dt
    v_new = v + ((AARE_dt./AARE_mass) .* neb_f); %% Euler velocity integrator

    neb_it = neb_it + 1;
    [Ynew,ind,covariance] = GPR_predict_pinv(X,p_new,values,alpha,optimizedParams,K);
    energy = Ynew(:,1);
    force = Ynew(:,2:end);
    tan=improved_tan(p_new,energy);
    spring = spring_parallel(p_new,tan,sp_const);
    [g_over,f_over] = neb_force(p_new,tan,spring,force);
    maxforce(neb_it) = max(max(g_over)).*con_fact; % maximum NEB force in eV/A
    fprintf(fileID, '%d\n', maxforce(neb_it));
    normforce(neb_it) = norm(g_over.*con_fact);
    % Write initial value to file
    fprintf(fileID1, '%d\n', normforce(neb_it));


    ant= (angle(g_over,CG_f));

    if ((ant>120))
        k=k+1;
        v=v/2;
        AARE_dt = AARE_dt * AARE_fdec;

        %%Euler integrator x_i+1 = x_i + v_i dt + 1/2 a_i dt^2
        p_new = points + (AARE_dt * v) + ((AARE_dt^2)*(neb_f)./(2*AARE_mass));

        %%% v_i+1 = v_i(mod) + Force/m dt
        v_new = v + ((AARE_dt./AARE_mass) .* neb_f);

        
        neb_it = neb_it + 1;
        [Ynew,ind,covariance] = GPR_predict_pinv(X,p_new,values,alpha,optimizedParams,K);
        energy = Ynew(:,1);
        force = Ynew(:,2:end);
        tan=improved_tan(p_new,energy);
        spring = spring_parallel(p_new,tan,sp_const);
        [g_new,f_new] = neb_force(p_new,tan,spring,force);
        maxforce(neb_it) = max(max(g_new)).*con_fact; % maximum NEB force in eV/A
        fprintf(fileID, '%d\n', maxforce(neb_it));
        normforce(neb_it) = norm(g_new).*con_fact;
    	% Write initial value to file
        fprintf(fileID1, '%d\n', normforce(neb_it));

    else

        g_new=g_over;
        f_new=f_over;
    end


    points=p_new; %%% assign the current position
    v=v_new;


    i=i+1; % Next iteration

    grad2 = grad1;
    grad1 = g_new;
    f = f_new;
    neb_f = grad1;
    max(f);
    
end

Ynew = Ynew*10000; %in kJ
final_normforce = normforce(neb_it);
disp(['Final NEB force = ', num2str(final_normforce)])
if isnan(final_normforce) || final_normforce == Inf
    neb_flag = 0;
else
    neb_flag = 1;
end
disp(['Number of NEB evaluations = ', num2str(neb_it)])
disp('NEB optimisation completed')

fclose(fileID);
fclose(fileID1);
toc
end
