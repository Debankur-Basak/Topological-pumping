function Delta = Delta_j(z, j, Delta_bar, Z)
    Delta = (-1)^j * Delta_bar * cos(2*pi*z / Z);
end