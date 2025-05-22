%% Function contains all constants and fluid parameters

function [mu_c,rho_c,rho_d,nu_c,D,A,L] = fun_constants(t)

mu_c = 1.846*10^(-5);
rho_c = 1.2;
rho_d = 1000;
nu_c = mu_c/rho_c;

D0 = 10^(-4);
lambda = 7*10^(-10);
D = sqrt(D0^2-lambda*t);

A = 50; % to be adapted
L = 2;

end