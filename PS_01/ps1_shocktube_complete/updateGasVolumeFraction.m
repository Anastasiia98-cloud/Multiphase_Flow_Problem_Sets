%% Function calculates the gas volume fraction from given vapor mass fraction, mixture density

function [alphav] = updateGasVolumeFraction(Yv,rho,p,T)

[Cv,~,gammav,~,~,~,pinfv,~,~,~] = stiffenedGasParameters;

rhov = (p+pinfv)/(Cv*T*(gammav-1)); % Vapor density

alphav = Yv*rho/rhov; % Calculation of gas volume fraction

end