% Script computes droplet diameter and mass over time for specified conditions
%% Specified conditions
Sh = 2;
Pr = 0.76;
Sc = 0.65;
lambda_c = 0.0257;
rho_d = 1000;
cp_c = 1005;
%% Plug in Yv_s and Yv_inf from problem 2
%Your Code Here
Yv_s = 0.00847;
Yv_inf = 0.00556;

%% Set initial Diameter
D_0 = 10^(-6);

%% Calculate lambda for Dsquared law
%Your Code Here
lambda = 4*Sh*Pr*lambda_c*(Yv_s-Yv_inf)/(Sc*rho_d*cp_c);

%% Calculate time limits
t_min = 0;
% Calculate final time
%Your Code Here
t_evap = D_0^2/lambda;
t = t_min:(t_evap-t_min)/1000:t_evap;

%% Dsquared law and mass calculation
%Your Code Here
D = sqrt(D_0^2-lambda*t);
m = (rho_d*pi*(D_0^2-lambda*t).^(3/2))/6;

%% Plots
figure(1);
subplot(2,1,1);
plot(t,D);
xlabel('t [s]');
ylabel('D [m]');

subplot(2,1,2);
plot(t,D.*D);
xlabel('t [s]');
ylabel('D^2 [m^2]');

figure(2);
plot(t,m);
xlabel('t [s]');
ylabel('m [kg]');


