clc, clear;

files = dir('*PROCESSED.mat');

load('./Geometry/Lp_m.mat');
load('./Geometry/Vf_m3.mat');

rho = 1.2;
c   = 340;

figure(1); clf;
axh(1) = nexttile(); hold on;% 
axh(3) = nexttile(); hold on;% 
axh(5) = nexttile(); hold on;% 

axh(2) = nexttile(); hold on;% 
axh(4) = nexttile(); hold on;% 
axh(6) = nexttile(); hold on;% 
for idx = 1 : length(files)
    
    filename = files(idx).name;
    load(filename);
    numpipe = str2num(filename(2:3));
    
    Vf = Vf_m3(numpipe);
    Lp = Lp_m(numpipe);
    
    scatter(axh(1), Area1*Vf, PRTfoot.*f1,'b');
    axh(1).XScale='log';axh(1).YScale='log';
    xlabel(axh(1),'PRTnondim');
    ylabel(axh(1),'Area1*Vf');
    
    scatter(axh(2), Area2*Vf, PRTfoot.*f1,'b');
    axh(2).XScale='log';axh(2).YScale='log';
    xlabel(axh(2),'PRTnondim');
    ylabel(axh(2),'Area2*Vf');
    
%     scatter(axh(3), Area1*Vf, betafit,'b');  % betafit/f1 ?
    scatter(axh(3), Area1*Vf, betafit.*PRTfoot,'b');  % betafit/f1 ?
    axh(3).XScale='log';axh(3).YScale='log';
    xlabel(axh(3),'beta');
    ylabel(axh(3),'Area1*Vf');
    
%     scatter(axh(4), Area2*Vf, betafit,'b'); % betafit/f1 ?
    scatter(axh(4), Area2*Vf, betafit.*PRTfoot,'b'); % betafit/f1 ?
    axh(4).XScale='log';axh(4).YScale='log';
    xlabel(axh(4),'beta');
    ylabel(axh(4),'Area2*Vf');
    
    scatter(axh(5), Area1*Vf, nufit,'b');
    axh(5).XScale='log';axh(5).YScale='log';
    xlabel(axh(5),'nu');
    ylabel(axh(5),'Area1*Vf');
    
    scatter(axh(6), Area2*Vf, nufit,'b');
    axh(6).XScale='log';axh(6).YScale='log';
    xlabel(axh(6),'nu');
    ylabel(axh(6),'Area2*Vf');
end

grid(axh, 'on');
box(axh, 'on');