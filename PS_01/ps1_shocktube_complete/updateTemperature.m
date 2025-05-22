%% Function computes Temperature from SG-EOS

function [T] = updateTemperature(rho,p)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,~,~,pinfv,pinfl,~,~] = stiffenedGasParameters;

%% Compute Temperature from SG-EOS
% Your Code Here
T = (p+pinfv)/((gammav-1)*Cv*rho);

end