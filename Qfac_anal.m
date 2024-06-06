clc, clear;

c = 340;

load('./Geometry/Lp_m.mat');
load('./Geometry/DpManipFitted_mm.mat');
load('./Geometry/f1manipfitted.mat');
load('./Geometry/f1manip.mat');

Dp_m = DpManipFitted_mm*1e-3;
rads = 0.5*Dp_m;
f1   = f1manipfitted;
f1   = f1manip;



mask = [3,5,10,13,15,17,24,25,29,34,37,39,41,44,46,48,51,53]';
Dp_m = Dp_m(mask);
Lp_m = Lp_m(mask);
rads = rads(mask);

n = 1;

chi   = 3e-5 * sqrt(n*f1)./rads;
OM    = n*pi*c./(Lp_m + 1.2*rads);
part1 = OM/(2*c);
part2 = (Lp_m + 1.2*rads)./(chi.*Lp_m + (OM.^2.*rads.^2)/(2*c^2) )  ;
Q1 = part1.*part2;

figure(3);clf;
plot(12*log2(f1/440), Q1,'o');
xlabel('$12 log_2 (f_1 / f_{ref})$', 'interpreter','latex');
ylabel('$Q-factor$', 'interpreter','latex');

