%%  Monte-Carlo sampling of particle injection
%% Constants
u0 = 100;
n_sample = 10000;

% Log-normal distribution of D
% Plug in given conditions
m =
nu =
mu =
sigma =
Dmin =
Dmax =
amin =
amax =

%% Sampling
%generate random numbers (column matrix, n_sample*1)
D = 
D = D(D>=Dmin & D<=Dmax); %Removes out of range values
count = length(D);
%generate random numbers (column matrix, count*1)
a_int =
a0 = (amax-amin)*a_int+amin;

%Empty arrays
D_start = D;
a0_start = a0;
[t_end,u_end,v_end,y_end] = deal(zeros(count,1));

% Calculate trajectory
for i=1:count
    [t,y_out] = particle_trajectory(u0,a0(i),D(i));
    t_end(i) = t(end);
    u_end(i) = y_out(end,1);
    v_end(i) = y_out(end,2);
    y_end(i) = y_out(end,4);
end

%% Plotting
figure(1);
subplot(1,2,1);
x = logspace(-4,-2,40);
histogram(D_start,x,'Normalization','probability');
xlabel('Size(m)'); ylabel('Count fraction'); title('Particle size distribution');
set(gca,'xscale','log');

subplot(1,2,2);
histogram(a0_start,'Normalization','probability');
xlabel('Angle'); ylabel('Count fraction'); title('Injection angle distribution');

u_t_end = [u_end t_end];
y_t_end = [y_end t_end];
y_u_end = [y_end u_end];
D_a0_start = [D_start a0_start];

figure(2);
subplot(2,2,1);
hist3(u_t_end);
xlabel('U velocity (m/s)'); ylabel('time (s)'); title('U velocity distribution');

subplot(2,2,2);
hist3(y_t_end);
xlabel('Y position (m)'); ylabel('time (s)'); title('Y position distribution');

subplot(2,2,3);
hist3(y_u_end);
xlabel('Y position (m)'); ylabel('U velocity (m/s)'); title('Y position with U velocity distribution');

subplot(2,2,4);
hist3(D_a0_start);
xlabel('Injection diameter (m)'); ylabel('Injection Angle'); title('Injection distributon');