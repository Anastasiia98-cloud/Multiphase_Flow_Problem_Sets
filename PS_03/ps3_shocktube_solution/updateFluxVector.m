%% Updates flux vector from primitive variables

function [F] = updateFluxVector(rho,Yv,u,p,E,Nx)
%% UPDATE ENTRIES OF FLUX VECTOR FROM PRIMITIVE VARIABLES
% Your Code Here
F(1,:) = rho(:,1).*u(:,1);
F(2,:) = rho(:,1).*u(:,1).*Yv(:,1);
F(3,:) = rho(:,1).*u(:,1).^2 + p(:,1);
F(4,:) = (E(:,1) + p(:,1)).*u(:,1);
end