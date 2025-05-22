%% Function updates internal energy

function [e] = updateIntEnergy(T,Yv,Yl,p)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

% Calculation of internal energy
e = Yl*(Cl*T*(p+gammal*pinfl)/(p+pinfl)+ql)+Yv*(Cv*T*(p+gammav*pinfv)/(p+pinfv)+qv);

end