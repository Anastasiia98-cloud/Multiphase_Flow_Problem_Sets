%% Script computes and plots trajectory of water droplet in air flow vortex flow field

% Flow field vectors
nx = 100; % to be adapted
ny = 100; % to be adapted
Lx = 10;
Ly = 10;
x = zeros(nx+1,1);
y = zeros(ny+1,1);
for i=1:(nx+1)
    x(i) = Lx/nx*i;
end
for i=1:(ny+1)
    y(i) = Ly/ny*i;
end

t_0 = 0; % to be adapted
t_f = 0.36;

ud_0 = 0;
vd_0 = 0;
xd_0 = 6;
yd_0 = 4;

u_c = zeros(nx+1,ny+1);
v_c = zeros(nx+1,ny+1);

% Constants / fluid properties
[mu_c,rho_c,rho_d,nu_c,~,~,~] = fun_constants();

for i=1:(nx+1)
    for j=1:(ny+1)
        [u_c(i,j),v_c(i,j)] = velocity_cont(x(j),y(i));
    end
end

quiver(x,y,u_c,v_c);

tspan = [t_0 t_f];
y0 = [ud_0 vd_0 xd_0 yd_0];
[t,y_out] = ode45(@odefun,tspan,y0);
hold on;
xlabel('x [m]');
ylabel('y x[m]');
plot(y_out(:,3),y_out(:,4));



