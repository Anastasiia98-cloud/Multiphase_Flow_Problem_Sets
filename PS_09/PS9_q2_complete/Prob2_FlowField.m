%% Script computes and plots trajectory of water droplet in air flow (Taylor-Green vortex)

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
t = t_0;
t_max = 5;
nt = 100;
dt = (t_max-t_0)/nt;

uc = zeros(nx+1,ny+1);
vc = zeros(nx+1,ny+1);

A = 50; % to be adapted
L = 2;

% Constants / fluid properties
mu_c = 1.846*10^(-5);
rho_c = 1.2;
nu_c = mu_c/rho_c;

for i=1:(nx+1)
    for j=1:(ny+1)
        uc(i,j) = A*cos(x(j)/L)*sin(y(i)/L);
        vc(i,j) = -A*sin(x(j)/L)*cos(y(i)/L);
    end
end

quiver(x,y,uc,vc);
xlabel('x [m]');
ylabel('y [m]');

