clc;
clear all;
close all;

%% System parameters
tau = 0.11;                   % in cm^-1
del_tilde = 1.5;              % dimensionless quantity
Delta_bar = 0.23;             % in cm^-1
Z = 74;                       % in cm 
pbc = 1;                      % Periodic Boundary Condition
N = 30;                       % No. of unit cells.
dz = 0.01;                    % Step size for z
z_values = (0:dz:Z);          % Range of z values

%% Initialize storage for eigenvalues
eigenvalues = zeros(2*N, length(z_values));

%% Eigenvalues and Eigenvectors at z = 0
[V_rs, D_rs] = eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc));
lambda = diag(D_rs);
[EPSA, ind] = sort(lambda);
V = V_rs(:, ind);

%% Projection Operator
P = zeros(2*N);
desired_band_start = 1;
desired_band_end = N;

for ii = desired_band_start:desired_band_end
    p = V(:, ii) * V(:, ii)';
    P = double(p + P);
end

D = zeros(2*N, 1);
for ii = 1:2*N
    D(ii, 1) = exp(1i * ii * 2 * pi / (2 * N));
end

X = diag(D);
Xp = P * X * P;

%% Wannier Function
[W, EV] = eig(Xp);

%% Disorder parameters
Xi_off = 0;
g_values = 0.1:0.1:0.4;  % g values to iterate over
Xi_on_values = 0:0.1:0.5; % Xi_on values to iterate over

%% Initialize storage for CM_shift
CM_shift = zeros(length(g_values), length(Xi_on_values));

%% Loop over Xi_on and g values
for jj = 1:length(g_values)
    g = g_values(jj);
    
    for kk = 1:length(Xi_on_values)
        Xi_on = Xi_on_values(kk);
        
        %% Finding the soliton
        itr = 0;         % Iterations
        max_itr = 501;   % Maximum number of iterations
        E_old = W(:, N/2); % Initial guess for the soliton state
        
        while itr <= max_itr
            L_H = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc);
            NL_H = g * diag(abs(E_old).^2);
            H_new = L_H - NL_H;
            
            [eig_vec, eig_val] = eig(H_new);
            lam = diag(eig_val);                     % Extract eigenvalues
            [E, ind] = sort(lam);                    % Sort eigenvalues
            V1 = eig_vec(:, ind);                    % Reorder eigenvectors
            
            band_indices = (desired_band_start:desired_band_end);
            eig_vec_band = V1(:, band_indices);
            eig_val_band = E(band_indices);
            
            [min_E, idx] = min(eig_val_band);       % Find the smallest eigenvalue
            E_new = eig_vec_band(:, idx);           % Extract the corresponding eigenvector
            
            RelativeTol = abs(sum(abs(E_new)) - sum(abs(E_old))) / sum(abs(E_new));
            
            if log10(RelativeTol) < -12
                break;
            else
                E_old = E_new;
            end
            
            itr = itr + 1;
        end
        
        %% Dynamics
        Ei = E_old;
        Reltolerance = 1e-12;
        N_periods = 2;
        PSI_sol = zeros(2*N, 1);
        
        % Load or generate disorder matrix
        % H_dis = Xi_on * eye(2*N); % Example: Identity matrix as disorder
        data = load('disordered_matrix.mat', 'H0');
        H_dis = Xi_on * data.H0;
        
        for Tt = 1:N_periods
            Tspan = [0, dz];
            option = odeset('RelTol', Reltolerance);
            
            for ii = 1:length(z_values) - 1
                H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc) + H_dis;
                
                [tt, E] = ode45(@(tt, E) Propagator(g, tt, E, H0), Tspan, Ei, option);
                E = E.';
                Ei = E(:, end);
                
                if ii == 1 && Tt == 1
                    PSI_sol = [PSI_sol, E(:, 1), Ei];
                else
                    PSI_sol = [PSI_sol, Ei];
                end
            end
        end
        
        PSI_sol = PSI_sol(:, 2:end);
        
        %% Center of Mass Calculation
        M = (1:2*N)';
        CM = M' * (abs(PSI_sol).^2);
        CM_shift(jj, kk) = abs(CM(end) - CM(1));
    end
end

%% Plot CM_shift vs Xi_on for different g values
figure;
hold on;
colors = lines(length(g_values)); % Generate distinct colors for each g value

for jj = 1:length(g_values)
    plot(Xi_on_values, CM_shift(jj, :)/N_periods,'o-', 'Color', colors(jj, :), 'LineWidth', 2, 'DisplayName', ['g = ', num2str(g_values(jj))]);
end



xlabel('Xi_{on}/ \tau');
xlim([0, 10]); % Set x-axis limits
ylabel('CM_{shift}');
title('CM_{shift} vs Xi_{on} for different g values');
legend show;
grid on;
hold off;