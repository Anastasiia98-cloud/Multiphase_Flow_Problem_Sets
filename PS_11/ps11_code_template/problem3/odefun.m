%% Function for ode45 in Prob3_Trajectory.m

function dydt = odefun(t,y,D)

[mu_c,~,rho_d,~,~] = fun_constants();

u_d = y(1);
v_d = y(2);
x_d = y(3);
y_d = y(4);

[u_c,v_c] = velocity_cont();

f1 = ffactor(u_c,v_c,u_d,v_d,D);

tauv = rho_d*D^2/(18*mu_c);

dydt = zeros(4,1);

dydt(1) = f1/tauv*(u_c-u_d);
dydt(2) = f1/tauv*(v_c-v_d);
dydt(3) = u_d;
dydt(4) = v_d;

end

function [u_c,v_c] = velocity_cont()

u_c=0;
v_c=10;

end

function f1 = ffactor(u_c,v_c,u_d,v_d,D)

[mu_c,rho_c,~,~,~] = fun_constants();

Re = rho_c*sqrt((u_c-u_d)^2+(v_c-v_d)^2)*D/mu_c;

if Re<1 && Re>=0
    f1 = 1;
elseif Re>=1 && Re<800
    f1 = 1+0.15*Re^0.687;
elseif Re>=800 && Re<2*10^5
    f1 = 1+0.15*Re^0.687+0.0175*Re/(1+42500*Re^(-1.16));
else
    disp('Re is out of bounds for f-factor');
end


end