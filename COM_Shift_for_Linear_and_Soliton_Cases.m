clc;
clear all;
close all;

%% System Parameters
tau = 0.55;                   % in cm^-1
del_tilde = 2.11;             % Dimensionless quantity
Delta_bar = tau * 6.97;       % in cm^-1
Z_max = 30 * 1 / tau;         % in cm 
pbc = 0;                      % Periodic Boundary Condition
N = 30;                       % Number of unit cells
N_periods = 3;                % Number of periods for propagation
dz = 0.01;                    % Step size for z
z_values = (0:dz:Z_max)*N_periods;        % Range of z values

Reltolerance = 1e-10;         % Relative tolerance for ODE solver


%% Wannier Function Initialization
[V_rs, D_rs] = eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z_max, 0, pbc));
lambda = diag(D_rs);
[~, ind] = sort(lambda);
V = V_rs(:, ind);

% Projection operator
P = zeros(2 * N);
desired_band_start = 1;
desired_band_end = N;
for ii = desired_band_start:desired_band_end
    p = V(:, ii) * V(:, ii)';
    P = double(p + P);
end

% Wannier function calculation
D = exp(1i * (1:2 * N) * 2 * pi / (2 * N)).';
X = diag(D);
Xp = P * X * P;
[W, ~] = eig(Xp);
W1 = W(:, 27); % Wannier function for site 27

%% Soliton Calculation

itr = 0;                     % Iterations
max_itr = 501;               % Maximum number of iterations
E_old = W1;                  % Initial guess for soliton
g = 1;                    % Nonlinearity parameter

while itr <= max_itr
    L_H = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z_max, 0, pbc);
    NL_H = g * diag(abs(E_old).^2); % Nonlinear Hamiltonian
    H_new = L_H - NL_H;             % Total Hamiltonian

    % Find eigenvalues and eigenvectors
    [eig_vec, eig_val] = eig(H_new);
    lam = diag(eig_val);
    [E, ind] = sort(lam);
    V1 = eig_vec(:, ind);

    % Extract the desired band
    band_indices = (desired_band_start:desired_band_end);
    eig_vec_band = V1(:, band_indices);
    eig_val_band = E(band_indices);

    % Find the smallest eigenvalue and its eigenvector
    [~, idx] = min(eig_val_band);
    E_new = eig_vec_band(:, idx);

    % Check relative tolerance
    RelativeTol = abs(sum(abs(E_new)) - sum(abs(E_old))) / sum(abs(E_new));
    if log10(RelativeTol) < -10
        break
    else
        E_old = E_new;
    end

    itr = itr + 1;
end


%% Initialize Storage for COM Shift Calculation.

Xi_on_range = linspace(0, 10, 40);                     % Xi_on values from 0 to 4
com_shift_linear = zeros(1, length(Xi_on_range));     % COM shift for linear case
com_shift_soliton = zeros(1, length(Xi_on_range));    % COM shift for soliton case

%% Loop over Xi_on values.


for Xi_on_idx = 1:length(Xi_on_range)
    Xi_on = Xi_on_range(Xi_on_idx);
    Xi_off = 0; % Set Xi_off = 0 for now

    %% Linear Case: Ei = W1
    M = (1:2 * N).'; % Position indices
    PSI_sol_linear = zeros(2 * N, length(z_values)); % State storage
    PSI_sol_linear(:, 1) = W1; % Initial state (linear)

    Ei = W1; % Initial state
    for ii = 1:length(z_values) - 1
        z = z_values(ii);
        H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z_max, z, pbc) + ...
             DisorderRiceMeleHamiltonian(N, pbc, Xi_on, Xi_off);
        [~, E] = ode45(@(tt, E) Linear_Propagator(tt, E, H0), [0 dz], Ei, odeset('RelTol', Reltolerance));
        Ei = E(end, :).'; % Update state for next step
        PSI_sol_linear(:, ii + 1) = Ei; % Store propagated state
    end

    % Compute center of mass for linear case
    normalized_PSI_sol_linear = abs(PSI_sol_linear).^2;
    CM_linear = M' * normalized_PSI_sol_linear;
    com_shift_linear(Xi_on_idx) = abs(CM_linear(end) - CM_linear(1));
    

    H_dis=DisorderRiceMeleHamiltonian(N, 0, Xi_on, Xi_off);

    %% Soliton Case: Ei = E_old
    PSI_sol_soliton = zeros(2 * N, length(z_values)); % State storage
    PSI_sol_soliton(:, 1) = E_old; % Initial state (soliton)

    Ei = E_old; % Initial state
    for ii = 1:length(z_values) - 1
        z = z_values(ii);
        H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z_max, z, pbc) + ...
             H_dis;
        [~, E] = ode45(@(tt, E) Propagator(g,tt, E, H0), [0 dz], Ei, odeset('RelTol', Reltolerance));
        Ei = E(end, :).'; % Update state for next step
        PSI_sol_soliton(:, ii + 1) = Ei; % Store propagated state
    end

    % Compute center of mass for soliton case
    normalized_PSI_sol_soliton = abs(PSI_sol_soliton).^2;
    CM_soliton = M' * normalized_PSI_sol_soliton;
    com_shift_soliton(Xi_on_idx) = abs(CM_soliton(end) - CM_soliton(1));
end

%% Plot COM Shift for Linear and Soliton Cases
figure;
plot(Xi_on_range, com_shift_linear, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Linear Case');
hold on;
plot(Xi_on_range, com_shift_soliton, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Soliton Case');
grid on;
xlabel('\xi_{on}');
ylabel('Center of Mass Shift');
title('COM Shift: Linear vs Soliton');
legend('show');
hold off;