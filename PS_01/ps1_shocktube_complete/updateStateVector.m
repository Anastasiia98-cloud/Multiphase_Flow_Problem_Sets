%% Updates state vector from primitive variables

function [W] = updateStateVector(rho,u,p,E,Nx)
%% UPDATE ENTRIES OF FLUX VECTOR FROM PRIMITIVE VARIABLES
% Your Code Here
W(1,:) = rho(:,1);
W(2,:) = rho(:,1) .* u(:,1);
W(3,:) = E(:,1);
end