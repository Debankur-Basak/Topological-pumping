function H0 = DisorderRiceMeleHamiltonian(N, pbc, Xi_on, Xi_off)
    % Initialize the Hamiltonian
    H0 = zeros(2*N);

    % Generate random disorder terms
    xi_on = 2 * rand(1, 2*N) - 1;    % On-site disorder [-1, 1]
    xi_off = 2 * rand(1, 2*N - 1) - 1; % Off-site disorder [-1, 1]

    for j = 1:2*N
        

        % On-site term including disorder
        H0(j, j) = Xi_on * xi_on(j);

        % Hopping terms including off-site disorder
        if j < 2*N
            hopping_term =Xi_off * xi_off(j);
            H0(j, j+1) = hopping_term;
            H0(j+1, j) = hopping_term;
        end
    end

    % Periodic boundary condition
    if pbc == 1
        
        hopping_pbc =Xi_off * xi_off(end);
        H0(1, 2*N) = hopping_pbc;
        H0(2*N, 1) = hopping_pbc;
    end
end