 H0 = zeros(2*N);

    
xi_on = 2 * rand(1, 2*N) - 1;    % On-site disorder [-1, 1]
xi_off = 2 * rand(1, 2*N - 1) - 1; % Off-site disorder [-1, 1]

for j = 1:2*N
    

    % On-site term including disorder
    H0(j, j) =xi_on(j);

    
end

save('disordered_matrix.mat', 'H0');
