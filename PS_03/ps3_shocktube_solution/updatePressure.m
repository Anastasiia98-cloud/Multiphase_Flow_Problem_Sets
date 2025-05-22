%% Function computes Pressure from SG-EOS

function [p] = updatePressure(rho,Yv,e)

% Definition of parameters from SG parameter vector
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,~,~] = stiffenedGasParameters;

%% Compute Pressure from SG-EOS
% Your Code Here
A = Yv*Cv+(1-Yv)*Cl;
B = rho.*(Yv*(gammav-1)*Cv+(1-Yv)*(gammal-1)*Cl).*(-e+Yv*qv+(1-Yv)*ql)+Yv*Cv*(pinfl+gammav*pinfv)+(1-Yv)*Cl*(pinfv+gammal*pinfl);
C = rho.*(Yv*(gammav-1)*Cv*pinfl+(1-Yv)*(gammal-1)*Cl*pinfv).*(-e+Yv*qv+(1-Yv)*ql)+Yv*Cv*gammav*pinfv*pinfl+(1-Yv)*Cl*gammal*pinfv*pinfl;
p = (-B+sqrt(B.^2-4*A.*C))./(2*A);

end