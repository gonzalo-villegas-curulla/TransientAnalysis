fs = 10e3;
dt = 1/fs;

T = 1;
t = [0:dt:T];
to = 0.4;

beta = 100;
nu = 0.01;

fun = 1./(1 + nu*exp(-beta*(t-to))).^(1/nu);
figure(1); clf;
plot(t, fun);

