%% Function calculates saturation temperature for given equilibrium pressure p
% Present implementation only works for two-phase mixture where both
% components in the phases have the pressure p

function [r_eq] = eqRadius(p,T,epsilong)

r_eq = 10^(-5); % Initial guess
[pv,pl] = getpressures(p,T,r_eq);
[delg] = gibbsEnergiesDirect(pv,pl,T);

if delg<0
    delta_r_eq = 10^(-7);
elseif delg>0
    delta_r_eq = -10^(-7);
end

if abs(delg)>epsilong
    
    delg_m1 =delg; % Definition of deltag from previous step (n=0)
    
    r_eq_m1 = r_eq;
    r_eq = r_eq+delta_r_eq; % Calculate r at step n=1
    
    [pv,pl] = getpressures(p,T,r_eq);
    [delg] = gibbsEnergiesDirect(pv,pl,T);
    
end

while abs(delg)>epsilong
    
    while delg*delg_m1<0 % Decrease iteration step 
        delta_r_eq = delta_r_eq/2;
        r_eq = r_eq_m1+delta_r_eq;
        [pv,pl] = getpressures(p,T,r_eq);
        [delg] = gibbsEnergiesDirect(pv,pl,T);
    end
    
    r_eq_m1 = r_eq;
    r_eq = r_eq+delta_r_eq; % Calculate Yv at next iterative step
    
    delg_m1 = delg; % Assign deltag from previous iterative step to deltag_m1
    [pv,pl] = getpressures(p,T,r_eq);
    [delg] = gibbsEnergiesDirect(pv,pl,T);

end

end