%% Function computes and plots trajectory of particle in air flow field

function [t,y_out] =particle_trajectory(u0,a0,D)

t_0 = 0; % to be adapted
t_f = 10;

ud_0 = 
vd_0 = 
xd_0 = 0;
yd_0 = 0;

% Constants / fluid properties
%[mu_c,rho_c,rho_d,nu_c,~] = fun_constants();

tspan = [t_0 t_f];
y0 = [ud_0 vd_0 xd_0 yd_0];

options = odeset('Events',@eventsFcn);
[t,y_out] = ode45(@(t,y) odefun(t,y,D),tspan,y0,options);


end


function [position,isterminal,direction] = eventsFcn(t,y)
position = y(3)-1; % Approaches 0 at x=L
isterminal = 1;
direction = 0;
end



