%% Rosin-Rammler computations

Dmin = 0;
Dmax = Inf;
delta = 144;
n = 2;

mu = integral(@(D)fun(D,delta,n),Dmin,Dmax);

var = integral(@(D)fun2(D,delta,n),Dmin,Dmax);
sigma = sqrt(var);

D = 1:1:1000;

for i=1:length(D)
    fm(i) = funPDF(D(i),delta,n);
    Fm(i) = funCDF(D(i),delta,n);
end

figure(1);
plot(D,fm);

figure(2);
plot(D,Fm);

%% Expressions
% Plug in expressions
function [fm] = funPDF(D,delta,n)
fm =
end

function [Fm] = funCDF(D,delta,n)
Fm =
end

function [fmD] = fun(D,delta,n)
fmD =
end

function [fmDsq] = fun2(D,delta,n)
fmDsq =
end