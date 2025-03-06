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

dz = 0.01;                    % Step size for z
z_values = 0:dz:Z_max;        % Range of z values

Reltolerance = 1e-10;         % Relative tolerance for ODE solver
N_periods = 1;                % Number of periods for propagation

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

%% Initialize Storage for COM Shift Calculation
Z_range = linspace(0, 40 * 1 / tau, 100); % Z values for COM shift
com_shift_values = zeros(1, length(Z_range)); % COM shift values

%% Compute Center of Mass Shift for Different Z
M = (1:2 * N).'; % Position indices

for z_idx = 1:length(Z_range)
    current_Z = Z_range(z_idx);
    z_values = 0:dz:current_Z; % Update z range

    % Initialize state for propagation
    PSI_sol = zeros(2 * N, length(z_values)); % State storage
    PSI_sol(:, 1) = W1; % Initial state

    % Propagation over z
    Ei = W1; % Initial state
    for ii = 1:length(z_values) - 1
        z = z_values(ii);
        H0 = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, current_Z, z, pbc);

        % Propagate using ODE solver
        [~, E] = ode45(@(tt, E) Propagator(tt, E, H0), [0 dz], Ei, odeset('RelTol', Reltolerance));
        Ei = E(end, :).'; % Update state for next step
        PSI_sol(:, ii + 1) = Ei; % Store propagated state
    end

    % Compute center of mass
    normalized_PSI_sol = abs(PSI_sol).^2;
    CM = M' * normalized_PSI_sol; % Compute center of mass at each z
    com_shift_values(z_idx) = abs(CM(end) - CM(1)); % COM shift for current Z
end

%% Plot COM Shift as a Function of Z
figure;
plot(Z_range, com_shift_values, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Z (cm)');
ylabel('Center of Mass Shift');
title('Center of Mass Shift as a Function of Z');