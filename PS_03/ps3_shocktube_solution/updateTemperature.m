%% Function computes Temperature from SG-EOS

function [T] = updateTemperature(rho,Yv,Yl,p)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,~,~,pinfv,pinfl,~,~] = stiffenedGasParameters;

%% Compute Temperature from SG-EOS
% Your Code Here
T = 1/(rho.*((gammal-1)*Yl*Cl./(p+pinfl)+(gammav-1)*Yv*Cv./(p+pinfv)));

end