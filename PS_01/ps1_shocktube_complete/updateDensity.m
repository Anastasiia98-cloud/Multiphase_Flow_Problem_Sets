function [rho] = updateDensity(p,T,Yv,Yl)

[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

% Calculation of gas and liquid densities
rhov = (p+pinfv)/((gammav-1)*Cv*T);
rhol = (p+pinfl)/((gammal-1)*Cl*T);

rho = 1/(Yv/rhov+Yl/rhol);


end