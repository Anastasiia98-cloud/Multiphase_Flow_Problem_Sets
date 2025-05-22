%% Function contains all constants and fluid parameters

function [mu_c,rho_c,rho_d,nu_c,Lx] = fun_constants()

mu_c = 1.846*10^(-5);
rho_c = 1.2;
rho_d = 1000;
nu_c = mu_c/rho_c;
Lx = 1;

end