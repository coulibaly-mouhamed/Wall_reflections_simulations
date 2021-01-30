function [rhs1, rhs2] = compute_rhs_full(phi_hat, eta_hat, t, p)
% Computes the RHS of the hyperbolic part
L =p.g(t);
for i=1:length(p.g(t))
    rhs1      =+   -L(i).*eta_hat+ p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat; 
end
rhs2      =   DtN(phi_hat,p) + 2*p.nu0.*p.K2.*eta_hat; % DtN computes phi_z_hat 
%rhs1 = rhs1 + p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat;
   
end