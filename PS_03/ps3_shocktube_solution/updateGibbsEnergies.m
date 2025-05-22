%% Function calculates the Gibbs energy differences of the liquid and vapor phases

function [deltag,meang,p,T] = updateGibbsEnergies(rho,e,Yv)

%% COMPUTE MIXTURE PRESSURE p, MIXTURE TEMPERATURE T, GIBBS FREE ENERGY DIFFERENCE deltag, MEAN VALUE meang
% Your Code Here
[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,qvp,qlp] = stiffenedGasParameters;

p = updatePressure(rho,Yv,e); 
Yl = 1-Yv;
T = updateTemperature(rho,Yv,Yl,p);

% Calculate Gibbs energies
gv = (gammav*Cv-qvp)*T-Cv*T*log(T^gammav/(p+pinfv)^(gammav-1))+qv;
gl = (gammal*Cl-qlp)*T-Cl*T*log(T^gammal/(p+pinfl)^(gammal-1))+ql;

deltag = gv-gl;
meang = (abs(gv)+abs(gl))/2;


%% %%

end