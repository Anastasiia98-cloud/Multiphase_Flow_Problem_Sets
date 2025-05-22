%% Function relaxes Gibbs free energy by means of iteration

function [Yv] = gibbsRelaxation(rho,e,Yv)

epsilon = 10^(-6); % Tolerance value that shall remain after total evaporation/condensation
epsilon_bar = 10^(-5); % Threshold value for activation of relaxation procedure 
epsilonv = 10^(-6); % Terminating parameter for bisection procedure
epsilon0 = 0.5; % Parameter for initialization of deltaYg
DY = 2; % Parameter for reinitialization of deltaYg in case of sign change of deltag

if Yv>epsilon_bar && Yv<1-epsilon_bar

    %% UPDATE GIBBS FREE ENERGY DIFFERENCE deltag AND MEAN VALUE meang
    
    
    
    %% %%

if deltag<0 % Evaporation takes place
    
    Yv_max = 1-epsilon;
    
    %% UPDATE MIXTURE PRESSURE p_max AND GIBBS FREE ENERGY DIFFERENCE deltag_max at Yv_max %%
    
    
    
    %% %%
    
    if deltag_max>0 || p_max<=0 % The actual volume fraction lies between the initial one and the maximum one (to exclude cases with negative pressure from direct allocation) // else: the maximum Yv is taken with the respective p and T
        Yv = bisectionGibbsEnergies(rho,e,Yv,deltag,meang,epsilon,epsilonv,epsilon0,DY,1);   % sign (last variable) indicates evap./cond.: for evaporation, sign=1
    else
        Yv = Yv_max; % Take maximum value of Yv
    end
        
elseif deltag>0 % Condensation takes place
    
    Yv_min = epsilon;
    
    %% UPDATE MIXTURE PRESSURE p_min AND GIBBS FREE ENERGY DIFFERENCE deltag_min at Yv_min %%
    
    
    
    %% %%
    
    if deltag_min<0 || p_min<=0 % The actual volume fraction lies between the minimum one and the initial one (to exclude cases with negative pressure from direct allocation) // else: the minimum Yv is taken with the respective p and T
        Yv = bisectionGibbsEnergies(rho,e,Yv,deltag,meang,epsilon,epsilonv,epsilon0,DY,-1);   % sign (last variable indicates evap./cond.: for condensation, sign=-1
    else
        Yv = Yv_min; % Take minimum value of Yv
    end
    
end

end

end