%% Updates primitive flow variables from state vector

function [rho,Yv,Yl,u,E,e,p,T,c] = updateFlowVariables(W,Nx)
    %% UPDATE PRIMITIVE FLOW VARIABLES FROM STATE VECTOR %% 
    % Your Code Here
    rho(1:Nx-1,1) = W(1,:);
    Yv(1:Nx-1,1) = W(2,:)./W(1,:);
    Yl(1:Nx-1,1) = ones(Nx-1,1)-Yv(1:Nx-1,1);
    u(1:Nx-1,1) = W(3,:)./W(1,:);
    E(1:Nx-1,1) = W(4,:);
    e           = E./rho-0.5*u.^2;

    %%    %%

    p = zeros(Nx-1,1);
    T = zeros(Nx-1,1);
    c = zeros(Nx-1,1);

    for j=1:Nx-1

        % Calculate mixture pressure
        p(j,1) = updatePressure(rho(j,1),Yv(j,1),e(j,1));

        % Update temperature
        T(j,1) = updateTemperature(rho(j,1),Yv(j,1),Yl(j,1),p(j,1));

        % Update speed of sound
        c(j,1) = updateSoundSpeed(p(j,1),rho(j,1),T(j,1),Yv(j,1),Yl(j,1));

    end

end