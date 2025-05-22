%% Function calculates saturation pressure for given equilibrium temperature T
% Present implementation only works for two-phase mixture where both
% components in the phases have the pressure p

function [psat] = saturationPressure(T,epsilong)

psat = 10^3; % Initial guess
[delg] = gibbsEnergiesDirect(psat,psat,T);

if delg<0
    delta_psat = 100;
elseif delg>0
    delta_psat = -100;
end

if abs(delg)>epsilong
    
    delg_m1 =delg; % Definition of deltag from previous step (n=0)
    
    psat_m1 = psat;
    psat = psat+delta_psat; % Calculate Tsat at step n=1
    
    [delg] = gibbsEnergiesDirect(psat,psat,T);
    
end

while abs(delg)>epsilong
    
    while delg*delg_m1<0 % Decrease iteration step 
        delta_psat = delta_psat/2;
        psat = psat_m1+delta_psat;
        [delg] = gibbsEnergiesDirect(psat,psat,T);
    end
    
    psat_m1 = psat;
    psat = psat+delta_psat; % Calculate Yv at next iterative step
    
    delg_m1 = delg; % Assign deltag from previous iterative step to deltag_m1
    [delg] = gibbsEnergiesDirect(psat,psat,T);

end

end