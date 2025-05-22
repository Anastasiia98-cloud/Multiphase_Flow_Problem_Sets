%% Function returns the continuous flow field velocity at a given location


function [u_c,v_c] = velocity_cont(x,y)

[~,~,~,~,~,A,L] = fun_constants();

u_c = A*cos(x/L)*sin(y/L);
v_c = -A*sin(x/L)*cos(y/L);


end