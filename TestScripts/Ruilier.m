% Ruilier

lam = 0.980;
w = 5.8/2;
D = 5.8;
f = 24.5;
Beta = pi/2*D/f*w/lam

a = 0.11;
b = 1.10;
rho = 2*((exp(-b^2)*(1-exp((b^2)*(1-a^2))))/(b*sqrt(1-a^2)))^2

f


f = w/(0.71*lam/D)