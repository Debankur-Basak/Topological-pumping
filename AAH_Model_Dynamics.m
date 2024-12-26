

clc;
clear all;
close all;

%% System Specifications.

pbc=1;                   % System with periodic boundary conditions.
delta=0;                 % On-site energy is 0.
N=100;                   % Total number of cells.

L=8000;                  % Total length of the sample in mm.
omega=2*pi/L;            % Frequency.
alpha=-(2*pi)/12;        % Initial Phase.

d=22/1000;               % Average separation between two waveguides in mm.
del=2/1000;              % The spatial modulation strength in mm.
g=0.2;

Np=300+1;                % No. of spatial points along the longitudinal direction.
Delz=L/(Np-1);           % Smallest part of length along the longitudinal direction.

z=0:Delz:L;              % Array of all the z-points.

x1=d+del*cos((2*pi/3)+omega*z+alpha);
x2=2*d+del*cos((2*pi/3)*2+omega*z+alpha);
x3=3*d+del*cos((2*pi/3)*3+omega*z+alpha);

distance_ab=zeros(length(z),1);
distance_bc=zeros(length(z),1);
distance_ca=zeros(length(z),1);

for ii=1:1:length(z)

    distance_ab(ii)=sqrt((x1(ii)-x2(ii))^2);
    distance_bc(ii)=sqrt((x2(ii)-x3(ii))^2);
    distance_ca(ii)=sqrt((-3*d+x3(ii)-x1(ii))^2);

end

%% Coupling constants.


Jab=27.21*exp(-330*distance_ab);
Jbc=27.21*exp(-330*distance_bc);
Jca=27.21*exp(-330*distance_ca);

% figure(8);
% plot(omega*z,Jab);
% hold on;
% plot(omega*z,Jbc);
% plot(omega*z,Jca);
% hold off;
% ylabel('Couplings (in mm^{-1})')
% xlabel('\Omega z');
% legend({'J_{ab}','J_{bc}','J_{ca}'});


%% Energy Bands.

% for ii=1:1:Np-1
% 
%     H0=AAH_Model_Hamiltonian(N,Jab(ii),Jbc(ii),Jca(ii),delta,pbc);
%     En=eig(H0);
% 
%     figure(1);
%     plot(omega*ii*Delz*ones(size(En,1),1),En,'b.');
%     xlabel('\Omega z');
%     ylabel('Energy');
%     title('Band Structure');
%     hold on;
% 
% end

%% Eigenvalues and Eigenvectors of the initial Hamiltonian.

[V_rs,D_rs]=eig(AAH_Model_Hamiltonian(N,Jab(1),Jbc(1),Jca(1),delta,pbc));
lambda=diag(D_rs);
[EPSA,ind]=sort(lambda);
V=V_rs(:,ind);

% figure(2);
% plot(1:1:3*N,EPSA,'r.');
% title('Sorted Energy');
% xlabel('Eigenvalue number');
% ylabel('Energy');

%% Projection Operator.

P=zeros(3*N);

for ii=1:1:N
    p=V(:,ii)*V(:,ii)';
    P=double(p+P);
end

D=zeros(3*N,1);

for ii=1:1:3*N
    D(ii,1)=exp(1i*ii*2*pi/(3*N));
end

X=diag(D);
Xp=P*X*P;

%% Wannier Function.

[W,EV]=eig(Xp);

%% Calculation of IPR.

% for ii=1:1:3*N
%     W1=W(:,ii);
%     IPR=sum(abs(W1).^4)/(sum(abs(W1).^2))^2;
%     figure(3);
%     plot(ii,IPR,'bo');
%     xlabel('site');
%     ylabel('IPR');
%     title('Inverse Participation Ratio');
%     hold on;
% end

%% Overlap

% fig4=figure(Name="Wannier functions.");
% figure(fig4);
% for ii = 1:1:3*N
%     W1 = W(:,ii);
% 
%     bar(W1.^2); % Bar plot of squared values
%     title(['Iteration: ', num2str(ii)]); % Display iteration number
%     xlabel('Index'); % Label for the x-axis
%     ylabel('W1^2'); % Label for the y-axis
%     pause(0.5); % Pause for 0.5 seconds to visualize each bar plot
% end
W1 = W(:,N/2);
ylabel('|W|^2');
xlabel('site no.');
title('Wannier function.');

Overlap=W1'*V/sqrt(sum(abs(W1).^2));
Overlap=abs(Overlap).^2;
Overlap=Overlap';

Overlap_N=zeros(N, 1); 
Overlap_N(1)=Overlap(1);

for jn=2:2:N-1
    Overlap_N(jn)=(Overlap(jn)+Overlap(jn+1))/2;
    Overlap_N(jn+1)=Overlap_N(jn);
end

Overlap_N(N)=Overlap(N); 
Overlap(1:N)=Overlap_N;

Overlap_Mean=mean(Overlap_N);

figure(4); yyaxis left
plot(EPSA, (1:1:3*N), 'bo'); 
ylabel('Eigenvalue number');
xlabel('Energy'); 
hold on;

figure(4); yyaxis right
plot(EPSA, Overlap, 'r.');
ylabel('Overlap');


%% Finding out the soliton.

itr=0;         % Iterations.
max_itr=501;   % Maximum number of iterations.
E_old=W1;
NL_sign=1;

fig2=figure(Name='Iterations',Position=[643.4,50.6,874.4,327.2]);


figure(fig2);
subplot(2,3,[1,4]);
bar(1:1:3*N,abs(E_old).^2/(E_old'*E_old));
title('Guessed state');
ylabel('|\psi_{guess}|^2');
xlabel('Sites');
set(gcf,'color','w');
fontname(fig2,"Arial");
fontsize(fig2,11,"points");


while itr<=max_itr

    L_H=AAH_Model_Hamiltonian(N,Jab(1),Jbc(1),Jca(1),delta,pbc);
    NL_H=g*diag(abs(E_old).^2);
    H_new=L_H-NL_H;
    [eig_vec,eig_val]=eig(H_new);
    [min_E,idx]=min(diag(eig_val));
    E_new=eig_vec(:,idx);
    RelativeTol=abs(sum(abs(E_new))-sum(abs(E_old)))/sum(abs(E_new));
    
    figure(fig2);
    subplot(2,3,[3 6]);
    bar(1:1:3*N,abs(E_new).^2/(E_new'*E_new));
    title(['After',num2str(itr),'iterations']);
    ylabel('|\psi_{soliton}|^2');
    xlabel('Sites');
    set(gcf,'color','w');
    fontname(fig2,"Arial");
    fontsize(fig2,11,"points");
    hold on;

    figure(fig2);
    subplot(2,3,2);
    plot(itr,min_E,'ro',MarkerSize=6);
    ylabel('Energy');
    xlabel('Iteration');
    set(gcf,'color','w');
    fontname(fig2,"Arial");
    fontsize(fig2,11,"points");
    hold on;

    figure(fig2);
    subplot(2,3,5);
    plot(itr,log10(RelativeTol),'bo',MarkerSize=6);
    ylabel('Tolerance');
    xlabel('Iteration');
    set(gcf,'color','w');
    fontname(fig2,"Arial");
    fontsize(fig2,11,"points");
    hold on;
    
    if log10(RelativeTol)<-15
        break
    else
        E_old=E_new;
    end

    itr=itr+1;
end

hold off;

%% Propagation of the soliton.
% 
% fig1=figure(Name='Propagation');
% z_max=100;
% z_pt = 20*round(z_max);
% Ei = E_old;
% 
% H0=AAH_Model_Hamiltonian(N,Jab(1),Jbc(1),Jca(1),delta,pbc);
% z_span=(linspace(0,z_max,z_pt))';
% opts = odeset('RelTol',1e-8);
% ODE1 = @(tt,E)(Propagator(g,tt,E,H0));
% [t,E] = ode45(ODE1,z_span,Ei,opts);
% I0 = abs(E).^2/(Ei'*Ei);
% 
% figure(fig1);
% imagesc([0 z_max],[1 3*N],I0');
% colormap hot;
% xlabel('z (mm)');
% ylabel('Sites');
% title('Soliton');
% set(gcf,'color','w');fontname(fig1,"Arial");fontsize(fig1,11,"points");



%% Dynamics.

Ei=E_old;
Reltolerance=10^-10;
N_periods=3;

PSI_sol=zeros(3*N,1);

for Tt=1:1:N_periods
    Zspan=[0,Delz];
    option=odeset('RelTol',Reltolerance);

    for t=1:1:Np-1
        H0=AAH_Model_Hamiltonian(N,Jab(t),Jbc(t),Jca(t),delta,pbc);
        [tt, E]=ode45(@(tt, E)Propagator(g,tt, E,  H0), Zspan, Ei, option); 
        E=E.';
        Ei=E(:,end);

        if t==1 && Tt==1
            PSI_sol=[PSI_sol,E(:,1),Ei];
        else
            PSI_sol=[PSI_sol,Ei];

        end


    end

end

PSI_sol=PSI_sol(:,2:end);

% PP=zeros(1,size(PSI_sol,2));
% for ii=1:1:size(PSI_sol,2)
%     PP(1,ii)=sum(abs(PSI_sol(:,ii)).^2);
% end
% figure(7);
% plot(PP);

%% Propagation.

fig5=figure(Name='Soliton Propagation',Position=[643.4,50.6,874.4,327.2]);

prop = (0:(Np-1)*N_periods) * Delz; % Correct
WG_no = 1:3*N;

% Plot the imagesc data
figure(fig5);
imagesc([prop(1), prop(end)], [WG_no(1), WG_no(end)], abs(PSI_sol));
colormap('hot');
ylabel('Waveguide number');
xlabel('Propagation distance (mm)');
title('Evolution of the initial state.');
% Add a secondary y-axis
ax1 = gca; % Get handle for the current axis
ax2 = axes('Position', ax1.Position, 'YAxisLocation', 'right', 'Color', 'none');

% Link x-axes and adjust ticks or labels
linkaxes([ax1, ax2], 'x'); % Link x-axes
ax2.XTick = []; % Disable x-ticks for the secondary axis
ax2.YDir = ax1.YDir; % Match y-direction
ax2.YLim = ax1.YLim; % Match y-limits

% Set custom Y-ticks with a gap of 3
step_size = 10;
custom_ticks = WG_no(1):step_size:WG_no(end);
ax2.YTick = custom_ticks; 
ax2.YTickLabel = custom_ticks; % Optional: Scale or modify labels

%% Center of Mass.

M=zeros(3*N, 1);
M(1:3*N, 1)=1:1:3*N;

CM=M'*(abs(PSI_sol).^2);


size(PSI_sol)

figure(11);
plot(prop, CM, 'ro'); ylabel('X_{cm}');
xlabel('Propagation distance (mm)'); 
title('Centre of mass.');


%% Non-Linear Schrodinger Equation.


function dEdz=Propagator(g,tt,E,H0)
dEdz=-1i*(H0-g*diag(abs(E).^2))*E;
end