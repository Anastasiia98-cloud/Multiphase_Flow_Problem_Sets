%% Compute vapour and liquid phase pressures due to surface tension
function [pv,pl] = getpressures(p,T,r)

k = 2.1*10^(-7); % Eoetvoes constant
Tc = 652;
M = 0.018015;
rho = 1000;
vm = M/rho;

% Surface tension
% Your Code Here
sigma = k*(Tc-T)/vm^(2/3);

% Phase pressure
% Your Code Here
pv = p;
pl = p+2*sigma/r;

end