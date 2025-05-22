%% Function calculates saturation pressure for given equilibrium temperature T
% Present implementation only works for two-phase mixture where both
% components in the phases have the pressure p

function [psat] = saturationPressure(T,epsilong)

psat = 10^3; % Initial guess
[delg] = gibbsEnergiesDirect(psat,T);

if delg<0
    delta_psat = 100;
elseif delg>0
    delta_psat = -100;
end

if abs(delg)>epsilong
    
    delg_m1 =delg; % Definition of deltag from previous step (n=0)
    
    psat_m1 = psat;
    psat = psat+delta_psat; % Calculate Tsat at step n=1
    
    [delg] = gibbsEnergiesDirect(psat,T);
    
end

while abs(delg)>epsilong
    
    while delg*delg_m1<0 % Decrease iteration step 
        delta_psat = delta_psat/2;
        psat = psat_m1+delta_psat;
        [delg] = gibbsEnergiesDirect(psat,T);
    end
    
    psat_m1 = psat;
    psat = psat+delta_psat; % Calculate Yv at next iterative step
    
    delg_m1 = delg; % Assign deltag from previous iterative step to deltag_m1
    [delg] = gibbsEnergiesDirect(psat,T);

end

end

function [delg] = gibbsEnergiesDirect(p,T)

[Cv,Cl,gammav,gammal,qv,ql,pinfv,pinfl,qvp,qlp] = stiffenedGasParameters;

% Calculate Gibbs energies
gv = (gammav*Cv-qvp)*T-Cv*T*log(T^gammav/(p+pinfv)^(gammav-1))+qv;
gl = (gammal*Cl-qlp)*T-Cl*T*log(T^gammal/(p+pinfl)^(gammal-1))+ql;

deltag = gv-gl;
meang = (abs(gv)+abs(gl))/2;
delg = deltag/meang;

end