%% Function computes Pressure from SG-EOS

function [p] = updatePressure(rho,e)

% Definition of parameters from SG parameter vector
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

%% Compute Pressure from SG-EOS
% Your Code Here
p = (gammav-1)*(e-qv)*rho - gammav*pinfv;

end