% Computes droplet evaporation time 
%% Specified Conditions
Sh = 2;
Nu = 2;
Pr = 0.76;
Sc = 0.65;
lambda_c = 0.0257;
rho_d = 1000;
cp_c = 1005;
L = 2454000;
xrel = 0.1;
M_H2O = 18;
M_M = 29;
epsilon = 10^(-1);

D_0 = 10^(-3);

T_min = 280;
T_max = 315;
p_min = 10^4;
p_max = 5*10^5;

%% Initialisation
T = T_min:1:T_max; nT=length(T);
p = p_min:10000:p_max; np = length(p);
t_evap = zeros(nT,np);

%% Calculate
%Your Code Here
for i=1:nT
    for j=1:np
        
        psat = saturationPressure(T(i),10^(-6));
        pref = psat;

        %Plug expression for Yv_inf at relative humidity level 
        Yv_inf = xrel*pref/p(j)*M_H2O/M_M;
        
        T_d = wetbulbTemperature(p(j),T(i),Sh,Nu,Pr,Sc,L,cp_c,Yv_inf,M_H2O,M_M,epsilon);
        
        %Plug expression for Yv_s at relative humidity level
        Yv_s = saturationPressure(T_d,10^(-6))/p(j)*M_H2O/M_M;
        
        %Plug expression to calculate lambda for D squared law
        lambda = 4*Sh*Pr*lambda_c*(Yv_s-Yv_inf)/(Sc*rho_d*cp_c);
        
        t_evap(i,j) = D_0^2/lambda;        
    end
end

%% Plotting
[X,Y]=meshgrid(p,T);
surf(X,Y,t_evap);
ylabel("Temperature (K)");
xlabel("Pressure (N/m^2)");
zlabel("t_{evap}")
