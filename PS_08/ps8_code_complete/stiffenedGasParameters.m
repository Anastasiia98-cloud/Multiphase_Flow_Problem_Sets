%% Parameters of liquid/vaporous water


function [Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,qvp,qlp] = stiffenedGasParameters
    
Cv = 1040;
Cl = 1816;
gammav = 1.43;
gammal = 2.35;
qv = 2030000;
ql = -1167000;
pinfv = 0;
pinfl = 10^9; 
qvp = -23.4*10^3;
qlp = 0;

end