%% Function solves the mixture pressure for an assumed Yv in the Gibbs energy relaxation process

function [p] = solvePQuadratic(Yv,rho,e)

% Definition of parameters from SG parameter vector
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

% Calculate parameters A, B, C of quadratic equation
A = Yv*Cv+(1-Yv)*Cl;
B = rho*(Yv*(gammav-1)*Cv+(1-Yv)*(gammal-1)*Cl)*(-e+Yv*qv+(1-Yv)*ql)+Yv*Cv*(pinfl+gammav*pinfv)+(1-Yv)*Cl*(pinfv+gammal*pinfl);
C = rho*(Yv*(gammav-1)*Cv*pinfl+(1-Yv)*(gammal-1)*Cl*pinfv)*(-e+Yv*qv+(1-Yv)*ql)+Yv*Cv*gammav*pinfv*pinfl+(1-Yv)*Cl*gammal*pinfv*pinfl;

%  Calculate the mixture pressure for the assumed Yv
p = (-B+sqrt(B^2-4*A*C))/(2*A);  % Only the positive sign of the sqrt expression is taken for resonable results

end
