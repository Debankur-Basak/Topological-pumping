
%% System parameters
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
N_periods = 2;
M = (1:2*N)';

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

D = exp(1i * (1:2*N) * 2 * pi / (2 * N)).';
X = diag(D);
Xp = P * X * P;

%% Wannier Function
[W, EV] = eig(Xp);


%% Disorder parameters
Xi_on = 0;
g_values = 0.1:0.1:0.7;  % g values to iterate over
Xi_off_values = 0:0.01:0.5; % Xi_on values to iterate over

%% Number of trials for averaging
num_trials = 10;

%% Initialize storage for CM_shift over trials
CM_shift_all = zeros(length(g_values), length(Xi_off_values), num_trials);
CM_shift_linear_all = zeros(length(Xi_off_values), num_trials);

%% Start parallel pool
if isempty(gcp('nocreate'))
    parpool; 
end

for trial = 1:num_trials
    %% Load or generate disorder matrix (new for each trial)
    %data = load('disordered_matrix.mat', 'H0');
   
    
    %% Parallel computation
    CM_shift = zeros(length(g_values), length(Xi_off_values));
    CM_shift_linear = zeros(1, length(Xi_off_values));

    parfor kk = 1:length(Xi_off_values)
        Xi_off = Xi_off_values(kk);
        H_dis = DisorderRiceMeleHamiltonian(N, 1, Xi_on, Xi_off);
        

        %% Dynamics for Linear Case
        Ei_linear = W(:, N/2);
        Reltolerance = 1e-12;
        PSI_sol_linear = zeros(2*N, 1);

        for Tt = 1:N_periods
            Tspan = [0, dz];
            option = odeset('RelTol', Reltolerance);

            for ii = 1:length(z_values) - 1
                H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc) + H_dis;
                [~, E_linear] = ode45(@(tt, E_linear) Propagator(0, tt, E_linear, H0), Tspan, Ei_linear, option);
                E_linear = E_linear.';
                Ei_linear = E_linear(:, end);

                if ii == 1 && Tt == 1
                    PSI_sol_linear = [PSI_sol_linear, E_linear(:, 1), Ei_linear];
                else
                    PSI_sol_linear = [PSI_sol_linear, Ei_linear];
                end
            end
        end

        PSI_sol_linear = PSI_sol_linear(:, 2:end);
        CM_linear = M' * (abs(PSI_sol_linear).^2);
        CM_shift_linear(kk) = abs(CM_linear(end) - CM_linear(1));

        %% Nonlinear Case
        CM_shift_local = zeros(length(g_values), 1);

        for jj = 1:length(g_values)
            g = g_values(jj);

            %% Finding the soliton
            itr = 0;
            max_itr = 501;
            E_old = W(:, N/2);

            while itr <= max_itr
                L_H = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc);
                NL_H = g * diag(abs(E_old).^2);
                H_new = L_H - NL_H;

                [eig_vec, eig_val] = eig(H_new);
                lam = diag(eig_val);
                [E, ind] = sort(lam);
                V1 = eig_vec(:, ind);

                band_indices = (1:N);
                eig_vec_band = V1(:, band_indices);
                eig_val_band = E(band_indices);

                [~, idx] = min(eig_val_band);
                E_new = eig_vec_band(:, idx);

                RelativeTol = abs(sum(abs(E_new)) - sum(abs(E_old))) / sum(abs(E_new));

                if log10(RelativeTol) < -12
                    break;
                else
                    E_old = E_new;
                end

                itr = itr + 1;
            end

            %% Dynamics for Soliton
            Ei = E_old;
            PSI_sol = zeros(2*N, 1);

            for Tt = 1:N_periods
                Tspan = [0, dz];
                option = odeset('RelTol', Reltolerance);

                for ii = 1:length(z_values) - 1
                    H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc) + H_dis;
                    [~, E] = ode45(@(tt, E) Propagator(g, tt, E, H0), Tspan, Ei, option);
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
            CM = M' * (abs(PSI_sol).^2);
            CM_shift_local(jj) = abs(CM(end) - CM(1));
        end

        CM_shift(:, kk) = CM_shift_local;
    end

    %% Store results for averaging
    CM_shift_all(:, :, trial) = CM_shift;
    CM_shift_linear_all(:, trial) = CM_shift_linear;
end

%% Close parallel pool
delete(gcp);

% %% Compute Mean CM_shift
% CM_shift_mean = mean(CM_shift_all, 3);
% CM_shift_linear_mean = mean(CM_shift_linear_all, 2);
% 
% %% Plot CM_shift vs Xi_on for different g values
% figure;
% hold on;
% colors = lines(length(g_values));
% 
% for jj = 1:length(g_values)
%     plot(Xi_off_values/tau, CM_shift_mean(jj, :)/N_periods, 'o-', 'Color', colors(jj, :), 'LineWidth', 2, ...
%          'DisplayName', ['g = ', num2str(g_values(jj))]);
% end
% 
% plot(Xi_off_values/tau, CM_shift_linear_mean/N_periods, 'o', 'DisplayName', 'Linear', 'Color', 'k');
% 
% xlabel('Xi_{on}/ \tau');
% ylabel('CM_{shift}');
% title('Averaged CM_{shift} vs Xi_{on} for different g values');
% legend show;
% grid on;
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute Mean and Standard Deviation of CM_shift
CM_shift_mean = mean(CM_shift_all, 3);
CM_shift_std = std(CM_shift_all, 0, 3);  % Standard deviation
CM_shift_linear_mean = mean(CM_shift_linear_all, 2);
CM_shift_linear_std = std(CM_shift_linear_all, 0, 2);

%% Plot CM_shift vs Xi_on with error bars
figure;
hold on;
colors = lines(length(g_values));

for jj = 1:length(g_values)
    % Error bar plot with '.' marking the mean
    errorbar(Xi_off_values/tau, CM_shift_mean(jj, :)/N_periods, CM_shift_std(jj, :)/N_periods, ...
        'o', 'Color', colors(jj, :), 'LineWidth', 1.5, 'MarkerFaceColor', colors(jj, :), ...
        'DisplayName', ['g = ', num2str(g_values(jj))]);
end

% Error bar for linear case
errorbar(Xi_off_values/tau, CM_shift_linear_mean/N_periods, CM_shift_linear_std/N_periods, ...
    'o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Linear');

xlabel('Xi_{on}/ \tau');
ylabel('CM_{shift}');
title('Averaged CM_{shift} with Error Bars vs Xi_{on}');
legend show;
grid on;
hold off;
