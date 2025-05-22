%% Function updates internal energy

function [e] = updateIntEnergy(T,p)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

% Calculation of internal energy
e = (p+gammav*pinfv)/(p+pinfv)*Cv*T+qv;

end