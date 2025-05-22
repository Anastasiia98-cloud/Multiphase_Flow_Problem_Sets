%% Function computes Speed of sound from SG-EOS

function [c] = updateSoundSpeed(p,rho,T,Yv,Yl)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,~,~,pinfv,pinfl,~,~] = stiffenedGasParameters;

%% Compute Speed of sound from SG-EOS
% Your Code Here
c = sqrt(gammav*(p+pinfv)/rho);

end