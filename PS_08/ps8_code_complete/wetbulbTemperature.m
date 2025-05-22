%% Function calculates wet-bulb temperature for given far-field temperature and pressure
% Present implementation only works for two-phase mixture where both
% components in the phases have the pressure p

function [T_d] = wetbulbTemperature(p,T_inf,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M,epsilon)

T_d = T_inf; % Initial guess
[delta] = solveEnergyEquation(T_d,T_inf,p,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M);

if delta<0
    delta_T = 10;
elseif delta>0
    delta_T = -10;
end

if abs(delta)>epsilon
    
    delta_m1 =delta; % Definition of delta from previous step (n=0)
    
    T_d_m1 = T_d;
    T_d = T_d+delta_T; % Calculate T_d at step n=1
    
    [delta] = solveEnergyEquation(T_d,T_inf,p,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M);
    
end

while abs(delta)>epsilon
    
    while delta*delta_m1<0 % Decrease iteration step 
        delta_T = delta_T/2;
        T_d = T_d_m1+delta_T;
        [delta] = solveEnergyEquation(T_d,T_inf,p,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M);
    end
    
    T_d_m1 = T_d;
    T_d = T_d+delta_T; % Calculate Td at next iterative step
    
    delta_m1 = delta; % Assign delta from previous iterative step to delta_m1
    [delta] = solveEnergyEquation(T_d,T_inf,p,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M);

end

end

function [delta] = solveEnergyEquation(T_d,T_inf,p,Sh,Nu,Pr,Sc,L,cp,Yv_c,M_H2O,M_M)

psat = saturationPressure(T_d,10^(-6));
if psat>p
    pref = p;
else
    pref = psat;
end

Yv_s = pref/p*M_H2O/M_M;

%Your Code Here
delta = T_d-T_inf+Sh/Nu*Pr/Sc*L/cp*(Yv_s-Yv_c);

end