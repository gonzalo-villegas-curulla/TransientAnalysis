clc, clear;

files = dir('*.mat');
nf = length(files);

mx = zeros(50, nf);

for idx = 1 : nf
    load(files(idx).name);
    mx(1:length(find(MX)),idx) = MX;
end


mx(mx==0) = nan;

figure(1); clf;
% boxplot(mx, 'PlotStyle', 'compact');
boxplot(1e3*mx, 'Widths',0.9);
ylabel('Downwards key moving time [ms]');
xlabel('Measured pipe number');