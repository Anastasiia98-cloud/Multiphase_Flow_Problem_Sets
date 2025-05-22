%% 1D Shock Tube Program
%
% W is the state vector
% F is the flux vector
% E/TE is the total energy
% SIE: Specific Internal Energy
% x : location of grid points
function [W,F,rho,u,p,T,E,e,x,time,frameAll] = shock_col

clc;
clear all;

% Length of domain (x-dir)
Lx = 1;

% Number of grid points
Nx = 1201;

% Final simulation time
final_time = 0.2;

% Cell size
hx = Lx/(Nx-1);

% Absolute location of cell centers
x = hx*[1:1:Nx-1]-hx/2;

% Initialize state and flux vectors to zero
W = zeros(3,Nx-1);
F = zeros(3,Nx-1);

% Initialize primitive flow variables and other parameters
rho = zeros(Nx-1,1); % Density
u = zeros(Nx-1,1);   % U-velocity
p = zeros(Nx-1,1);   % Pressure
T = zeros(Nx-1,1);   % Temperature
c = zeros(Nx-1,1);   % Mixture speed of sound
E = zeros(Nx-1,1);   % Total energy
e = zeros(Nx-1,1);   % Specific internal energy
CFL = 0.5;           % CFL number

% Initial conditions :1D shock tube

p(1:(Nx-1)/2,1) = 1;
p((Nx-1)/2:(Nx-1),1) = 0.1;

rho(1:(Nx-1)/2,1) = 1;
rho((Nx-1)/2:(Nx-1),1) = 0.125;

u(1:(Nx-1)/2,1) = 0;
u((Nx-1)/2:(Nx-1),1) = 0;

% Update temperature
for j=1:Nx-1
    T(j,1) = updateTemperature(rho(j,1),p(j,1));
end

% Update internal energy
for j=1:Nx-1
    e(j,1) = updateIntEnergy(T(j,1),p(j,1));
end

% Update mixture speed of sound
for j=1:Nx-1
    c(j,1) = updateSoundSpeed(p(j,1),rho(j,1));
end

E(1:Nx-1,1) = rho(1:Nx-1,1).*(e(1:Nx-1,1)+0.5*u(1:Nx-1,1).^2); % Update total energy

%% UPDATE STATE AND FLUX VECTORS %%
[W] = updateStateVector(rho,u,p,E,Nx);
[F] = updateFluxVector(rho,u,p,E,Nx);

%%   %%

% Compute time step
dt = CFL*hx/max(c + abs(u));

% Initial time, time step
time = 0;
ts = 1;

load('reference_data.mat');

% Loop over time until final time
while(time<final_time)
    
    % Solve for state vector
    [W] = solveEulerHLLC(W,F,p,c,Nx,hx,dt);
    
    %% MASS TRANSFER SUBSTEP AFTER HYPERBOLIC STEP -- in Problem Session 3! %%
    % MASS TRANSFER BY MEANS OF EQUILIBRIUM-BASED GIBBS FREE ENERGY RELAXATION 
    
    
    %% %%

    % Update primitive flow variables 
    [rho,u,E,e,p,T,c] = updateFlowVariables(W,Nx);
    
    % Update flux vector
    [F] = updateFluxVector(rho,u,p,E,Nx);
    
    % Advance time
    time = time + dt;
    
    % Update time step
    dt = CFL*hx/max(c + abs(u));
    
    % What follows is only for plotting / visualization
    
    display(['Time step: ' num2str(ts)]);
    display(['Current time : ' num2str(time)]);
    
%     if mod(ts,20) == 0
%         f1=figure(1);
%         semilogy(x,rho(:),'b','LineWidth',1);
%         hold on;
%         plot(x_ref,rho_ref,'r');
%         xlim([0 1]);
%         %ylim([1 10^3]);
%         xlabel('x');
%         ylabel('Density');
%         legend('Simulation','Reference');
%         hold off;
%         framerho(ts/20)=getframe(f1);
%     end
    
%     if mod(ts,20) == 0
%         f2=figure(2);
%         plot(x,Yv(:),'b','LineWidth',1);
%         hold on;
%         plot(x_ref,Yv_ref,'r');
%         xlim([0 1]);
%         ylim([0 1]);
%         xlabel('x');
%         ylabel('Vapor Mass Fraction');
%         legend('Simulation','Reference');
%         hold off;
%         frameY(ts/20)=getframe(f2);
%     end
    
%     if mod(ts,20) == 0
%         f3=figure(3);
%         semilogy(x,p(:),'b','LineWidth',1);
%         hold on;
%         plot(x_ref,p_ref,'r');
%         xlim([0 1]);
%         %ylim([10^5 10^8]);
%         xlabel('x');
%         ylabel('Pressure');
%         legend('Simulation','Reference');
%         hold off;
%         framep(ts/20)=getframe(f3);
%     end
    
%     if mod(ts,20) == 0
%         f4=figure(4);
%         plot(x,u(:),'b','LineWidth',1);
%         hold on;
%         plot(x_ref,u_ref,'r');
%         xlim([0 1]);
%         %ylim([0 350]);
%         xlabel('x');
%         ylabel('Velocity');
%         legend('Simulation','Reference');
%         hold off;
%         frameu(ts/20)=getframe(f4);
%     end
    
%     if mod(ts,20) == 0
%         f5=figure(5);
%         plot(x,T(:),'b','LineWidth',1);
%         hold on;
%         %plot(x_ref,T_ref,'r');
%         xlim([0 1]);
%         %ylim([500 1100]);~
%         xlabel('x');
%         ylabel('Temperature');
%         %legend('Simulation','Reference');
%         hold off;
%         frameT(ts/20)=getframe(f5);
%     end

    if mod(ts,20) == 0
        f6=figure(6);
        
        subplot(2,2,1);
        semilogy(x,rho(:),'b','LineWidth',1);
        hold on;
        plot(x_ref,rho_ref,'r');
        xlim([0 1]);
        %ylim([1 10^3]);
        xlabel('x');
        ylabel('Density');
        legend('Simulation','Reference');
        hold off;
        
%         subplot(2,3,2);
%         plot(x,Yv(:),'b','LineWidth',1);
%         hold on;
%         plot(x_ref,Yv_ref,'r');
%         xlim([0 1]);
%         ylim([0 1]);
%         xlabel('x');
%         ylabel('Vapor Mass Fraction');
%         legend('Simulation','Reference');
%         hold off;
        
        subplot(2,2,2);
        semilogy(x,p(:),'b','LineWidth',1);
        hold on;
        plot(x_ref,p_ref,'r');
        xlim([0 1]);
        %ylim([10^5 10^8]);
        xlabel('x');
        ylabel('Pressure');
        legend('Simulation','Reference');
        hold off;
        
        subplot(2,2,3);
        plot(x,u(:),'b','LineWidth',1);
        hold on;
        plot(x_ref,u_ref,'r');
        xlim([0 1]);
        %ylim([0 350]);
        xlabel('x');
        ylabel('Velocity');
        legend('Simulation','Reference');
        hold off;
        
        subplot(2,2,4);
        plot(x,T(:),'b','LineWidth',1);
        hold on;
        %plot(x_ref,T_ref,'r');
        xlim([0 1]);
        %ylim([500 1100]);~
        xlabel('x');
        ylabel('Temperature');
        %legend('Simulation','Reference');
        hold off;
        
        frameAll(ts/20)=getframe(f6);
    end
    
    ts = ts + 1;
    
end


end

