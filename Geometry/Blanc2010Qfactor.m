clc, clear;

xi  = 3e-5;
rho = 1.2;
co  = 340;

load('f1ManipMean.mat')
% load('f1manipfitted.mat')
load('Lp_m.mat');
load('DpManipFitted_mm.mat')

DpManip = 1e-3*DpManipFitted_mm;


% Blanc2010, Q2 =====================================
syms R n f1 L funpGVC(R,n) fun(R,n,f1,L)

fun(R,n,f1,L)  =  n*pi*0.5./(xi*sqrt(n*f1)*L/R + 0.5*(R*n*pi)^2./(L+1.2*R)^2);
funp = diff(fun, R);
funpGVC(R,n,f1,L) = -0.5*pi*n* (-xi*sqrt(n*f1)/R^2 + (n*R*pi)^2*( ((L+1.2*R)/R -1.2)/(L+1.2*R)^3 ) )./ (xi^2*(n*f1)/R^2  +...
    0.25*(n*pi*R)^4/(L+1.2*R)^4  + (xi*sqrt(n*f1)*R*(n*pi)^2)/(L+1.2*R)^2 );

Ls  = [0.04 : 0.005 : 1.5]';
f1s = co./(2*Ls);
Rs  = [0.001:0.0001:0.150]';

Rres1 = zeros(size(Ls));
Rres2 = zeros(size(Ls));
Rres3 = zeros(size(Ls));
% Rres4 = zeros(size(Ls));

fprintf('Progress: %1.2f',0.0);
    
for idx = 1 : length(Ls)
    
    thef1 = f1s(idx);
    theL  = Ls(idx);
    
%     jdx = find(funpGVC(Rs,2,thef1, theL)>0, 1,'last');  
    
    jdx = find(funp(Rs,1,thef1, theL)>0, 1,'last');    
    Rres1(idx) = Rs(jdx);

    jdx = find(funp(Rs,2,thef1, theL)>0, 1,'last');    
    Rres2(idx) = Rs(jdx);
    
    jdx = find(funp(Rs,3,thef1, theL)>0, 1,'last');    
    Rres3(idx) = Rs(jdx);
    
    %jdx = find(funp(Rs,4,thef1, theL)>0, 1,'last');    
    %Rres4(idx) = Rs(jdx);
    
    fprintf('\b\b\b\b%1.2f', idx/length(Ls));
end

Dres1 = 2*Rres1;
Dres2 = 2*Rres2;
Dres3 = 2*Rres3;
% Dres4 = 2*Rres4;

% Plot =====================================
%%
fo = 440;
LGD = {'Measured pipes','$Q_{max,1}$','$Q_{max,2}$','$Q_{max,3}$'};%,'$Q_4$'};
FSZ = 14; 
figure(1); clf; hold on; grid;
% mask = [3,5,10,13,15,17,24,25,29,34,37,39,41,44,44,46,48,51,53];
% mask  = [3,5,10,13,15,17,24,25,27, 29,32, 34,37,39,41,44]';
mask  = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27, 29,32, 34,37,39,41,44]';
mask  = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27, 29,32, 34,37,39,41,44,46,48,51,53]';
F1MEAN = [F1MEAN; 1760;1975;2349;2637];
plot( 12*log2(F1MEAN/fo), log(DpManip(mask)) ,'d');
%
plot( 12*log2(f1s/fo) , log(Dres1) ,'-k');
plot( 12*log2(f1s/fo) , log(Dres2) ,'-k');
plot( 12*log2(f1s/fo) , log(Dres3) ,'-k');
%plot( 12*log2(f1s/fo) , log10(Dres4) ,'-k');
xlabel('Frequency $12 log_2 (f_1/f_{440})$','interpreter','latex','fontsize',FSZ);
ylabel('$log$ (Diam [m])','interpreter','latex','fontsize',FSZ);
% legend(LGD,'location','northoutside','interpreter','latex');
legend(LGD,'location','best','interpreter','latex','fontsize',FSZ);
