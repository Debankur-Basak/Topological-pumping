function delta = del_j(z, j, tau, delta_tilde, Z)
    delta = tau * (exp((-1)^j * delta_tilde * sin(2*pi*z / Z)) - 1);
end

