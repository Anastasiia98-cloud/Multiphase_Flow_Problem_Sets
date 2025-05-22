%% Function to compute g difference between phases
function [delg] = gibbsEnergiesDirect(pv,pl,T)

[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,qvp,qlp] = stiffenedGasParameters;

% Calculate Gibbs energies
gv = (gammav*Cv-qvp)*T-Cv*T*log(T^gammav/(pv+pinfv)^(gammav-1))+qv;
gl = (gammal*Cl-qlp)*T-Cl*T*log(T^gammal/(pl+pinfl)^(gammal-1))+ql;

deltag = gv-gl;
meang = (abs(gv)+abs(gl))/2;
delg = deltag/meang;

end