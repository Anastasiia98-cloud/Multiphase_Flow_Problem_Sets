%% Script computes equilibrium radius for various pressures & temperatures
clc; clear;
T_min = 273;
T_max = 650;
p_min = 10000;
p_max = 10^7;
T = zeros(T_max-T_min+1,1);
p = zeros((p_max-p_min)/10000+1,1);
r = zeros(T_max-T_min+1,(p_max-p_min)/10000+1);
num_points = size(p,1)*size(T,1);
p_sat=zeros(size(T));
prog=0;

for i=T_min:T_max
    T(i-T_min+1,1) = i;
end
for i=p_min:10000:p_max
    p((i-p_min)/10000+1,1) = i;
end

epsilong = 10^(-6);
for i=1:length(T)
    p_sat(i,1) = saturationPressure(T(i,1),epsilong);
    if(mod(i,round(length(T)/10))==0)
        disp(['Calculating Saturation pressures: ' num2str(round(i/length(T)*100)) '% done']);
        pause(0.001);
    end
end

for i=T_min:T_max
   for j=p_min:10000:p_max
    prog=prog+1;
    if (p((j-p_min)/10000+1)>p_sat(i-T_min+1,1))
        r(i-T_min+1,(j-p_min)/10000+1) = eqRadius(p((j-p_min)/10000+1),T(i-T_min+1,1),epsilong);
    else
        r(i-T_min+1,(j-p_min)/10000+1) = NaN;
    end
    
    if(mod(prog,round(num_points/10))==0)
        disp(['Calculating radii: ' num2str(round(100*prog/num_points)) '% done']);
        pause(0.001);
    end
   end
end


%% Plotting
[X,Y] = meshgrid(p,T);
rnm = r*1e9;
h=gca;
surf(X,Y,rnm,'EdgeColor', 'none');
xlabel('Pressure (Pa)');
ylabel('Temperature (K)');
zlabel('Radius(nm)');
set(h,'zscale','log');
colormap('hot');
caxis([0.03 1]);