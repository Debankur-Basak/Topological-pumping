clc;
clear all;
close all;

%% System parametres

tau = 0.08;                   % in cm^-1
del_tilde = 2.11;             % dimentionless quantity
Delta_bar_values = linspace(0,tau*10,20);         % in cm^-1
CM_shift_values=zeros(size(Delta_bar_values));
for xx=1:1:length(Delta_bar_values)
    Delta_bar=Delta_bar_values(xx);
Z = 74;                    % in cm 
pbc=1;                        % Periodic Boundary Condition
N=30;                         % No. of unit cells.

dz = 0.01;                    % Step size for z

z_values = (0:dz:Z);          % Range of z values

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

site=15;


W1 = W(:,site); 


% figure(6);
% subplot(1,2,1);
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
% for ii=1:1:2*N
% 
% W1=W(:,ii);
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
% subplot(1,2,1);
% bar(W1.^2);
% ylabel('|W|^2');
% xlabel('site no.');
% title('Wannier function.');
% 
% subplot(1,2,2);
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
% pause(0.5);
% title("Iteration = ",num2str(ii));
% 
% end

%% Disorder parametres.

Xi_on=0;
Xi_off=0;
g=0;


%% Finding out the soliton.

% 
% itr=0;         % Iterations.
% max_itr=501;   % Maximum number of iterations.
% 
% E_old=W1;
% 
% 
% NL_sign=1;
% 
% fig8=figure(Name='Iterations',Position=[643.4,50.6,874.4,327.2]);
% 
% 
% figure(fig8);
% subplot(4,3,[1,4]);
% bar(1:1:2*N,abs(E_old).^2/(E_old'*E_old));
% title('Guessed state');
% ylabel('|\psi_{guess}|^2');
% xlabel('Sites');
% set(gcf,'color','w');
% fontname(fig8,"Arial");
% fontsize(fig8,11,"points");
% 
% subplot(4,3,[7,10]);
% site=27;
% site_range = site-5:site+5;
% angle_filtered = angle(W1(site_range));
% grid on;
% plot(site_range, angle_filtered / pi, 'o', 'LineWidth', 1.5);
% xlabel('Site no.');
% ylabel('Phase in radians/\pi');
% title('Phase of Wannier Function');
% grid off;
% 
% 
% while itr<=max_itr
% 
%     L_H=NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, 0, pbc);
%     %NL_H=DisorderRiceMeleHamiltonian(N, pbc, Xi_on, Xi_off);
%     NL_H=g*diag(abs(E_old).^2);
%     H_new = L_H - NL_H;
% 
% 
% 
%     [eig_vec, eig_val] = eig(H_new);
%     lam = diag(eig_val);                     % Extract eigenvalues
%     [E, ind] = sort(lam);                    % Sort eigenvalues
%     V1 = eig_vec(:, ind);                    % Reorder eigenvectors to match sorted eigenvalues
% 
%     band_indices = (desired_band_start:desired_band_end);
%     eig_vec_band = V1(:, band_indices);
%     eig_val_band =E(band_indices);
% 
% 
% 
% 
%     [min_E, idx] = min(eig_val_band);     % Find the smallest eigenvalue in the sorted list
%     E_new = eig_vec_band(:, idx);        % Extract the corresponding eigenvector
% 
% 
%     RelativeTol=abs(sum(abs(E_new))-sum(abs(E_old)))/sum(abs(E_new));
% 
%     figure(fig8);
%     subplot(4,3,[3 6]);
%     bar(1:1:2*N,abs(E_new).^2/(E_new'*E_new));
%     title(['After',num2str(itr),'iterations']);
%     ylabel('|\psi_{soliton}|^2');
%     xlabel('Sites');
%     set(gcf,'color','w');
%     fontname(fig8,"Arial");
%     fontsize(fig8,11,"points");
%     hold on;
% 
%     figure(fig8);
%     subplot(4,3,[9 12]);
%     site_range = site-5:site+5;
%     angle_filtered1 = angle(E_new(site_range));
%     grid on;
%     plot(site_range, angle_filtered1 / pi, 'o', 'LineWidth', 1.5);
%     xlabel('Site no.');
%     ylabel('Phase in radians/\pi');
%     title('Phase of Soliton');
%     grid off;
%     hold on;
% 
% 
% 
% 
% 
%     figure(fig8);
%     subplot(4,3,[2,5]);
%     plot(itr,min_E,'ro',MarkerSize=6);
%     ylabel('Energy');
%     xlabel('Iteration');
%     set(gcf,'color','w');
%     fontname(fig8,"Arial");
%     fontsize(fig8,11,"points");
%     hold on;
% 
%     figure(fig8);
%     subplot(4,3,[8,11]);
%     plot(itr,log10(RelativeTol),'bo',MarkerSize=6);
%     ylabel('Tolerance');
%     xlabel('Iteration');
%     set(gcf,'color','w');
%     fontname(fig8,"Arial");
%     fontsize(fig8,11,"points");
%     hold on;
% 
%     if log10(RelativeTol)<-12
%         break
%     else
%         E_old=E_new;
%     end
% 
%     itr=itr+1;
% end
% 
% hold off;


%% Dynamics.


%Ei=E_old;

Ei=W1;

% single_site=zeros(2*N,1);
% single_site(site,1)=1;
% Ei=single_site;

Reltolerance=10^-12;

N_periods=1;
PSI_sol=zeros(2*N,1);

%H_dis=DisorderRiceMeleHamiltonian(N, pbc, Xi_on, Xi_off);
%H_dis=Xi_on*eye(2*N);

data=load('disordered_matrix.mat', 'H0');
H_dis=data.H0;
H_dis=Xi_on*H_dis;

for Tt=1:1:N_periods
    Tspan=[0,dz];
    option=odeset('RelTol',Reltolerance);

    for ii=1:1:length(z_values)-1
        H0=NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z_values(ii), pbc);
            %+H_dis;
        
        [tt, E]=ode45(@(tt, E)Propagator(0,tt, E,  H0), Tspan, Ei, option); 
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

    








%% Initialize storage for overlap matrix


% OL = zeros(length(z_values), 2 * N); % Overlap storage

%% Loop over z values


% for zz = 1:length(z_values)
%     z = z_values(zz);
% 
%     % Eigenvectors for phi (periodic boundary conditions)
%     [phi_V_rs, phi_D_rs] = eig(NewRiceMeleHamiltonian(N, tau, del_tilde, Delta_bar, Z, z, 1));
%     phi_lambda = diag(phi_D_rs);
%     [~, phi_ind] = sort(phi_lambda);
%     phi_V = phi_V_rs(:, phi_ind); % Sorted eigenvectors
% 
%     psi = PSI_sol(:, zz); % Propagated state at current z
% 
%     % Compute overlap matrix for current z
%     for ii = 1:2*N
%         OL(zz, ii) = abs(psi' * phi_V(:, ii)); % Overlap magnitude
% 
%     end
% 
% 
% end

%% Plot the overlap matrix


% figure;
% imagesc(1:2*N, z_values/Z, OL'); % Transpose OL for correct alignment
% colorbar;
% ylabel(colorbar, 'Overlap Magnitude |< \Psi(z) | \Phi(z) >|');
% 
% % Set axis properties
% set(gca, 'YDir', 'normal'); % Ensure z=0 starts at the bottom
% xlabel('State Index, n');
% ylabel('Adiabatic Parameter z/Z');
% title('Overlap of Propagated State with Eigenstates');
% 
% % Adjust ticks for readability
% xticks(1:10:2*N); % Show every 10th state index
% yticks(0:0.2:1);  % Show z/Z in steps of 0.2


%% Propagation.

% fig9=figure(Name='Propagation',Position=[643.4,50.6,874.4,327.2]);
% 
% prop = (0:(length(z_values)-1)*N_periods) * dz; 
% WG_no = 1:2*N;
% 
% % Plot the imagesc data
% figure(fig9);
% imagesc([prop(1), prop(end)], [WG_no(1), WG_no(end)], abs(PSI_sol).^2);
% colormap('hot');
% ylabel('Waveguide number');
% xlabel('Propagation distance (mm)');
% title('Evolution of the initial state.');






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

%% Center of Mass.

M=zeros(2*N, 1);
M(1:2*N, 1)=1:1:2*N;

CM=M'*(abs(PSI_sol).^2);
CM_shift=abs(CM(end)-CM(1));
CM_shift_values(xx)=CM_shift;



end
figure(11);
grid on;
plot(Delta_bar_values/tau, CM_shift_values, 'ro-'); 
ylabel('X_{cm}');
xlabel('\Delta'); 
title('Centre of mass shift with \Delta.');
grid off;
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