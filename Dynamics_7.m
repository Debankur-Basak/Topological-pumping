clc;
clear all;
close all;

%% System parametres

tau = 0.55;                   % in cm^-1
del_tilde = 2.11;             % dimentionless quantity
Delta_bar = tau*6.97;         % in cm^-1
Z = 30*1/tau;                    % in cm 
pbc=1;                        % Periodic Boundary Condition
N=30;                         % No. of unit cells.

dz = 0.01;                    % Step size for z

z_values = (0:dz:Z);          % Range of z values

M(1:2*N, 1)=1:1:2*N;
%% Initialize storage for eigenvalues.

eigenvalues = zeros(2*N, length(z_values));

%% Calculate eigenvalues for each z.

% for idx = 1:length(z_values)
%     z = z_values(idx);
% 
%     H = NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z, pbc);
%     eigenvalues(:, idx) = eig(H); % Compute eigenvalues
% end

%% Plot eigenvalues as a function of z

% figure;
% hold on;
% for n = 1:2*N
%     plot(z_values/Z, eigenvalues(n, :), 'b.');
% end
% xlabel('z/Z');
% ylabel('Energy E (cm^{-1})');
% grid on;



%% Initialize storage for overlap matrix

% OL = zeros(length(z_values), 2 * N); % Overlap storage

%% Loop over z values.

% 
% for zz = 1:length(z_values)
%     z = z_values(zz);
% 
%     % Eigenvectors for psi (open boundary conditions)
%     [psi_V_rs, psi_D_rs] = eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z, 0));
%     psi_lambda = diag(psi_D_rs);
%     [~, psi_ind] = sort(psi_lambda);
%     psi_V = psi_V_rs(:, psi_ind); % Sorted eigenvectors
% 
%     % Eigenvectors for phi (periodic boundary conditions)
%     [phi_V_rs, phi_D_rs] = eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z, 1));
%     phi_lambda = diag(phi_D_rs);
%     [~, phi_ind] = sort(phi_lambda);
%     phi_V = phi_V_rs(:, phi_ind); % Sorted eigenvectors
% 
%     % Compute overlap matrix for current z
%     % OL(zz, ii) is the overlap between psi_V(:, ii) and phi_V(:, ii)
%     for ii = 1:2 * N
%         OL(zz, ii) = abs(psi_V(:, ii)' * phi_V(:, ii)); % Dot product (overlap)
%     end
% end



%% Plot the overlap matrix with proper z-axis orientation.


% figure;
% 
% % Transpose OL for proper alignment and plot
% imagesc(1:2*N, z_values / Z, OL'); 
% colorbar;
% 
% % Set axis labels and title
% xlabel('state number, n');
% ylabel('z / Z');
% title('< \Psi(z) | \Phi(z)>');
% 
% % Ensure z starts from bottom-left and increases upwards
% set(gca, 'YDir', 'normal'); % Normal direction for y-axis (bottom to top)
% 
% % Adjust tick marks (optional, for better readability)
% set(gca, 'YTick', linspace(0, max(z_values / Z), 5), ...
%          'YTickLabel', num2str(linspace(0, 1, 5)', '%.2f'));



%% Eigenvalues and Eigenvectors 

[V_rs,D_rs]=eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc));
lambda=diag(D_rs);
[EPSA,ind]=sort(lambda);
V=V_rs(:,ind);

% figure(3);
% plot(1:1:2*N,EPSA,'r.');
% title('Sorted Energy at time=0');
% xlabel('Eigenvalue number');
% ylabel('Energy');

%% Projection Operator.

P=zeros(2*N);
desired_band_start=1;
desired_band_end=N;

for ii=desired_band_start:desired_band_end
    p=V(:,ii)*V(:,ii)';
    P=double(p+P);
end

D=zeros(2*N,1);

for ii=1:1:2*N
    D(ii,1)=exp(1i*ii*2*pi/(2*N));
end

X=diag(D);
Xp=P*X*P;

%% Wannier Function.

[W,EV]=eig(Xp);

%% Calculation of IPR.

% for ii=1:1:2*N
%     W1=W(:,ii);
%     IPR=sum(abs(W1).^4)/(sum(abs(W1).^2))^2;
%     figure(4);
%     plot(ii,IPR,'bo');
%     xlabel('site');
%     ylabel('IPR');
%     title('Inverse Participation Ratio');
%     hold on;
% end

%% Checking for the wannier functions.


% for ii = 1:1:2*N
%     figure(5);
%     W1 = W(:,ii);
% 
%     bar(W1.^2); % Bar plot of squared values
%     title(['Iteration: ', num2str(ii)]); % Display iteration number
%     xlabel('Index'); % Label for the x-axis
%     ylabel('W1^2'); % Label for the y-axis
%     pause(1); % Pause for 0.5 seconds to visualize each bar plot
% end
% 
% 



%% Plotting the Intensity and phase of the wannier function.

site=27;
% figure(6);
% subplot(1,2,1);
W1 = W(:,site); 
% bar(W1.^2);
% ylabel('|W|^2');
% xlabel('site no.');
% title('Wannier function.');
% 
% subplot(1,2,2);
% site_range = site-5:site+5;
% angle_filtered = angle(W1(site_range)); 
% plot(site_range, angle_filtered / pi, 'o', 'LineWidth', 1.5);
% xlabel('Site no.');
% ylabel('Phase in radians/\pi');
% title('Phase of Wannier Function');
% grid on;

%% Overlap 
% 
% Overlap=W1'*V/sqrt(sum(abs(W1).^2));
% Overlap=abs(Overlap).^2;
% Overlap=Overlap';
% 
% Overlap_N=zeros(N, 1); 
% Overlap_N(1)=Overlap(1);
% 
% for jn=2:2:N-1
%     Overlap_N(jn)=(Overlap(jn)+Overlap(jn+1))/2;
%     Overlap_N(jn+1)=Overlap_N(jn);
% end
% 
% Overlap_N(N)=Overlap(N); 
% Overlap(1:N)=Overlap_N;
% 
% Overlap_Mean=mean(Overlap_N);
% 
% figure(7); 
% yyaxis left;
% plot(EPSA, (1:1:2*N), 'bo'); 
% ylabel('Eigenvalue number');
% xlabel('Energy'); 
% hold on;
% 
% figure(7);
% yyaxis right
% plot(EPSA, Overlap, 'ro');
% ylabel('Overlap');
% 


%% Disorder parametres.

Xi_on=2;
Xi_off=0;
g_values = [0.5, 1, 5, 10, 15, 20];
Xi_on_values = 0:1:15; 
CM_shift_soliton_g = zeros(length(g_values), length(Xi_on_values));

%CM_shift_linear = zeros(size(Xi_on_values));




for g_idx = 1:length(g_values)
    g = g_values(g_idx);

% Finding out the soliton.


itr=0;         % Iterations.
max_itr=501;   % Maximum number of iterations.

E_old=W1;


NL_sign=1;




while itr<=max_itr

    L_H=NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc);
   
    NL_H=g*diag(abs(E_old).^2);
    H_new = L_H - NL_H;



    [eig_vec, eig_val] = eig(H_new);
    lam = diag(eig_val);                     % Extract eigenvalues
    [E, ind] = sort(lam);                    % Sort eigenvalues
    V1 = eig_vec(:, ind);                    % Reorder eigenvectors to match sorted eigenvalues

    band_indices = (desired_band_start:desired_band_end);
    eig_vec_band = V1(:, band_indices);
    eig_val_band =E(band_indices);




    [min_E, idx] = min(eig_val_band);     % Find the smallest eigenvalue in the sorted list
    E_new = eig_vec_band(:, idx);        % Extract the corresponding eigenvector


    RelativeTol=abs(sum(abs(E_new))-sum(abs(E_old)))/sum(abs(E_new));

    

    if log10(RelativeTol)<-12
        break
    else
        E_old=E_new;
    end

    itr=itr+1;
end









for idx = 1:length(Xi_on_values)
    
    Xi_on = Xi_on_values(idx);



% Dynamics for the soliton.


Ei=E_old;

%Ei=W1;

% single_site=zeros(2*N,1);
% single_site(site,1)=1;
% Ei=single_site;

Reltolerance=10^-12;

N_periods=3;
PSI_sol=zeros(2*N,1);


data=load('disordered_matrix.mat', 'H0');
H_dis=data.H0;
H_dis=Xi_on*H_dis;




for Tt=1:1:N_periods
    Tspan=[0,dz];
    option=odeset('RelTol',Reltolerance);

    for ii=1:1:length(z_values)-1
        H0=NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc)+... 
            H_dis;
        %H0(145,145)=detune;
        [tt, E]=ode45(@(tt, E)Propagator(g,tt, E,  H0), Tspan, Ei, option); 
        E=E.';
        Ei=E(:,end);

        if ii==1 && Tt==1
            PSI_sol=[PSI_sol,E(:,1),Ei];
        else
            PSI_sol=[PSI_sol,Ei];

        end


    end

end

    PSI_sol=PSI_sol(:,2:end);

    
    CM_soliton = M' * (abs(PSI_sol).^2);
    CM_shift_soliton_g(g_idx, idx) = abs(CM_soliton(end) - CM_soliton(1)) / N_periods;


    % CM_shift = abs(CM(end) - CM(1)) / 3;
    % CM_final_values = [CM_final_values, CM_shift];


%end

%% Dynamics for the Linear case.

% Ei_lin=W1;
% PSI_sol_lin=zeros(2*N,1);

% for Tt=1:1:N_periods
%     Tspan=[0,dz];
%     option=odeset('RelTol',Reltolerance);
% 
%     for ii=1:1:length(z_values)-1
%         H0=NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc)+... 
%             H_dis;
%         %H0(145,145)=detune;
%         [tt, E]=ode45(@(tt, E)Propagator(0,tt, E,  H0), Tspan, Ei_lin, option); 
%         E=E.';
%         Ei_lin=E(:,end);
% 
%         if ii==1 && Tt==1
%             PSI_sol_lin=[PSI_sol_lin,E(:,1),Ei_lin];
%         else
%             PSI_sol_lin=[PSI_sol_lin,Ei_lin];
% 
%         end
% 
% 
%     end
% 
% end

    % PSI_sol_lin=PSI_sol_lin(:,2:end);
    % CM_linear = M' * (abs(PSI_sol_lin).^2);
    % CM_shift_linear(idx) = abs(CM_linear(end) - CM_linear(1)) / N_periods;

    end

end







%% Plot CM Shifts
figure;
hold on;

% Plot soliton CM shifts for each g value
colors = lines(length(g_values));
for g_idx = 1:length(g_values)
    plot(Xi_on_values, CM_shift_soliton_g(g_idx, :), 'o-', 'LineWidth', 1.5, ...
         'Color', colors(g_idx, :), 'DisplayName', sprintf('Soliton (g = %.1f)', g_values(g_idx)));
end

% % Plot linear case CM shift
% plot(Xi_on_values, CM_shift_linear, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Linear');

xlabel('\Xi_{on}');
ylabel('CM Shift');
title('Center of Mass Shift vs \Xi_{on}');
legend('Location', 'best');
grid on;
set(gcf, 'Color', 'w');























% figure;
% plot(Xi_on_values, CM_shift_soliton, 'ro', 'LineWidth', 1.5, 'DisplayName', 'Soliton');
% hold on;
% %plot(Xi_on_values, CM_shift_linear, 'bo', 'LineWidth', 1.5, 'DisplayName', 'Linear');
% xlabel('\Xi_{on}');
% ylabel('CM Shift');
% title('Center of Mass Shift vs \Xi_{on}');
% legend('Location', 'best');
% grid on;
% set(gcf, 'color', 'w');
% 




% Propagation 

% fig9=figure(Name='Propagation and CM',Position=[643.4,50.6,874.4,327.2]);
% 
% prop = (0:(length(z_values)-1)*N_periods) * dz; % Correct
% WG_no = 1:2*N;
% 
% % Plot the imagesc data
% figure(fig9);
% subplot(2,2,1);
% imagesc([prop(1), prop(end)], [WG_no(1), WG_no(end)], abs(PSI_sol).^2);
% colormap('hot');
% ylabel('Waveguide number');
% xlabel('Propagation distance (mm)');
% title('Evolution of the soliton');
% 
% 
% figure(fig9);
% subplot(2,2,2);
% 
% M=zeros(2*N, 1);
% M(1:2*N, 1)=1:1:2*N;
% grid on;
% CM=M'*(abs(PSI_sol).^2);
% 
% plot(prop, CM, 'ro'); ylabel('X_{cm}');
% xlabel('z'); 
% title('Centre of mass of the soliton.');
% grid off;
% 
% 
% 
% figure(fig9);
% subplot(2,2,3);
% imagesc([prop(1), prop(end)], [WG_no(1), WG_no(end)], abs(PSI_sol_lin).^2);
% colormap('hot');
% ylabel('Waveguide number');
% xlabel('Propagation distance (mm)');
% title('Evolution of the linear state');
% 
% 
% 
% 
% figure(fig9);
% subplot(2,2,4);
% grid on;
% CM_lin=M'*(abs(PSI_sol_lin).^2);
% 
% plot(prop, CM_lin, 'ro'); ylabel('X_{cm}');
% xlabel('z'); 
% title('Centre of mass of the linear state.');
% grid off;










% % Add a secondary y-axis
% ax1 = gca; % Get handle for the current axis
% ax2 = axes('Position', ax1.Position, 'YAxisLocation', 'right', 'Color', 'none');
% 
% % Link x-axes and adjust ticks or labels
% linkaxes([ax1, ax2], 'x'); % Link x-axes
% ax2.XTick = []; % Disable x-ticks for the secondary axis
% ax2.YDir = ax1.YDir; % Match y-direction
% ax2.YLim = ax1.YLim; % Match y-limits
% 
% % Set custom Y-ticks with a gap of 3
% step_size = 10;
% custom_ticks = WG_no(1):step_size:WG_no(end);
% ax2.YTick = custom_ticks; 
% ax2.YTickLabel = custom_ticks; % Optional: Scale or modify labels





%% Testing
% 
% for ii=1:100:size(PSI_sol,2)
%     figure;
%     hold on;
%     bar(abs(PSI_sol(:,ii)));
%     pause(0.5);
% 
% end
% 
% hold off;