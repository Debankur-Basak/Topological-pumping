function dEdz=Propagator(g,tt,E,H0)
dEdz=-1i*(H0-g*diag(abs(E).^2))*E;
end