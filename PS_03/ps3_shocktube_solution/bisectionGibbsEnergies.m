%% Function determines Yv with bisection method

function [Yv] = bisectionGibbsEnergies(rho,e,Yv,deltag,meang0,epsilon,epsilonv,epsilon0,DY,sign)

% Calculate initial deltaYv
if sign>0 % Evaporation
    deltaYv = epsilon0*(1-epsilon-Yv);
elseif sign<0 % Condensation
    deltaYv = epsilon0*(epsilon-Yv);
end

delg = deltag/meang0; % Calculate delg at n=1

if abs(delg)>epsilonv
    
    deltag_m1 =deltag; % Definition of deltag from previous step (n=0)
    Yv_m1 = Yv;
    Yv = Yv+deltaYv; % Calculate Yv at step n=1
    if Yv<epsilon  % To exclude cases where Yv gets out of bounds
        Yv = epsilon;
        deltaYv = Yv-Yv_m1;
    elseif Yv>1-epsilon
        Yv = 1-epsilon;
        deltaYv = Yv-Yv_m1;
    end
    
    p = updatePressure(rho,Yv,e);
    while p<=0  % To prevent that p gets negative with a too high deltaYv
        deltaYv = deltaYv/DY;
        Yv = Yv_m1+deltaYv;
        p = updatePressure(rho,Yv,e);
    end
    
    %% UPDATE GIBBS FREE ENERGY DIFFERENCE deltag %%
    % Your Code Here
    
    [deltag,~,~,~] = updateGibbsEnergies(rho,e,Yv);
    
    %% %%
    
    delg = deltag/meang0;
    
end

while abs(delg)>epsilonv
    
    % Define new deltaYv based on the behavior of deltag at current and previous iterative step (=bisection method)
    if deltag*deltag_m1<0
        deltaYv = -deltaYv/DY;
    elseif deltag*deltag_m1>0
        deltaYv = deltaYv/DY;
    end
    
    Yv_m1 = Yv;
    Yv = Yv+deltaYv; % Calculate Yv at next iterative step
    if Yv<epsilon  % To exclude cases where Yv gets out of bounds
        Yv = epsilon;
        deltaYv = Yv-Yv_m1;
    elseif Yv>1-epsilon
        Yv = 1-epsilon;
        deltaYv = Yv-Yv_m1;
    end
    
    p = updatePressure(rho,Yv,e);
    while p<=0  % To prevent that p gets negative with a too high deltaYv
        deltaYv = deltaYv/DY;
        Yv = Yv_m1+deltaYv;
        p = updatePressure(rho,Yv,e);
    end
    
    deltag_m1 = deltag; % Assign deltag from previous iterative step to deltag_m1
    
    %% UPDATE GIBBS FREE ENERGY DIFFERENCE deltag %%
    % Your Code Here
    
    [deltag,~,~,~] = updateGibbsEnergies(rho,e,Yv);
    
    %% %%
    
    delg = deltag/meang0; % Calculate delg at next iterative step
    %display(['deltag:' int2str(deltag)]);
end

end