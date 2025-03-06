function H0 = NewRiceMeleHamiltonian(N, tau, delta_tilde, Delta_bar, Z, z, pbc)

H0 = zeros(2*N);

for j = 1:2*N
    delta = del_j(z, j, tau, delta_tilde, Z);
    Delta = Delta_j(z, j, Delta_bar, Z);

    % Hopping terms
    if j < 2*N
        H0(j, j+1) = tau/2 + delta/2;
        H0(j+1, j) = tau/2 + delta/2;
    end

    % On-site terms
    H0(j, j) = Delta;
end

% Periodic boundary condition
if pbc==1
    H0(1, 2*N) = tau/2 + del_j(z, 1, tau, delta_tilde, Z)/2;
    H0(2*N, 1) = tau/2 + del_j(z, 2*N, tau, delta_tilde, Z)/2;
end

end