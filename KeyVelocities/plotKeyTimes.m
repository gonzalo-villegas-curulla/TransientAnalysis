clc, clear;
load('KoenigPalletWindWidth_m.mat');

files = dir('A*.mat');
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
ylim([0 35]);


figure(2); clf; hold on;
for idx = 1 : length(files)
    scatter(idx*0+1,1e3*mx(:,idx));
end
ylim([0 35]);
%%
pipenum = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27,29,32,34,37,39,41,44];
mxWithPalletWidth = mx;
for idx = 1 : length(files)
    pref = files(idx).name;
    pipe = str2num(pref(2:3));
    pallwid = palletwinwidth_m(idx);
    mxWithPalletWidth(:,idx) = (pallwid)*( 1./mxWithPalletWidth(:,idx)  );
end

figure(3); clf;
boxplot(mxWithPalletWidth, 'Widths',0.9);
ylabel('Downwards key moving time [ms]');
xlabel('Measured pipe number');
ylim([0 1.1]);