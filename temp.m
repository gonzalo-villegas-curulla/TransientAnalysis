clc, clear;

fs = 10e3;
dt = 1/fs;

T = 0.5;
t = [0:dt:T];
to = 0.25;

beta = 100;
nu   = 1e-7;


fun = 100./(1 + nu*exp(-beta*(t-to))).^(1./nu);

figure(12); clf;
plot(t, fun); grid on;