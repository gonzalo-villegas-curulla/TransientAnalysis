%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear;


% ===================================
% Options & Vena-Contracta (in/out)
% ===================================
VCin  = 0.95;  % Def. 0.36
VCout = 0.62;  % Def. 0.62
mouthbool = 0; % Take mouth measured pressure into account or not

% ===================================
%   Physical constants and params
% ===================================

rho = 1.2;              % [kg/m^3]
c   = 340; c2 = c^2;    % [m/s]
V4  = 0.084e-3;         % [m^3] (Def. 0.084e-3)
l4  = 1e-6;             % [m]  
l5  = 5e-6;             % [m] 
S4  = pi * (0.0020)^2;  % [m^2] measured ~0.0019m
S5  = 3.2940e-5;        % [m^2]

% ===================================
% Load test data 
% ===================================
load('Jussieu_Big_20221210_172043.mat');
clc;

fs     = 51.2e3; dt = 1/fs;

Tend   = PR_data.Time(end);
p3t    = PR_data.Time;              % Measured data time vector
p3     = PR_data.pressureData(:,1); % P_reservoir(t) measured
% p3     = p3-min(p3)*1.01;           % Avoid 
p3 = max(abs(p3),eps);


% win = hann(200);
% p3(1:100) = win(1:100).*p3(1:100);




[b,a] = butter(2, 0.25);
p3 = filtfilt(b,a,p3);

p4meas = PR_data.pressureData(:,2); % P_foot(t) measured
pmmeas = PR_data.pressureData(:,3);

% Time-Shift pipe pressure as if it were
% delayed mouth pressure plus an amplitude factor:
pmmeas = circshift(pmmeas,fix(fs*0.0005));

% ===================================
% DOWNSAMPLE DATA measured @51.2kHz
% ===================================
fsnew  = 4e3; dtnew = 1/fsnew; 
p3t    = resample(p3t,fsnew,fs); 
p3     = resample(p3,fsnew,fs);
p4meas = resample(p4meas,fsnew,fs);
pmmeas = resample(pmmeas,fsnew,fs);
fs     = fsnew; dt = dtnew;


% ===================================
%       SEARCH f1 and T1 
%    (in mouth pressure signal)
% ===================================

facoust = 340/(2*Lp);
try
    Ryin = yin(pmmeas,fs);
    f1   = Ryin.f0_Hz;
    N    = length(f1);
    ll   = fix(N/2)-fix(N*0.2) : fix(N/2)+fix(N*0.2);
    f1   = mean(f1(ll),'omitnan');
catch
    f1   = 285;
end

T1  = 1/f1;        % Fundamental period
Tn  = ceil(T1/dt); % Fund. period in samples



% ===================================
%      Re-adjust Sin and Sout
% ===================================

S4 = S4*VCin;
S5 = Hm*h*VCout;


% ===================================
%     Pre-compute coefficients
% ===================================

a1 = 1/(rho*l4);
a2 = -1/(rho*l4);
a3 = -1/(2*l4);
b1 = (S4*rho*c2)/V4;
b2 = -(S5*rho*c2)/V4;
c1 = 1/(rho*l5);
c2 = -1/(rho*l5);
c3 = -1/(2*l5);
params = [a1,a2,a3,b1,b2,c1,c2,c3];


% ===================================
% Prepare ODE integration conditions
% ===================================

tspan = [1e-4 Tend];
ic    = 1e-5*[1,1,1]; % <<<<< Put the initial values you have for y_k(0)

opts  = odeset('RelTol',1e-4,'AbsTol',1e-4,'OutputFcn',@odetpbar);
fprintf('Starting ODE solver...\n');
tic
[t,y] = ode15s(@(t,y) myodes(t,y,p3t,p3,pmmeas,params,mouthbool), tspan, ic, opts);
% [t,y] = ode113(@(t,y) myodes(t,y,p3t,p3,pmmeas,params), tspan, ic, []);
fprintf('(done)\n');
toc

% ===================================
%       Read-out solutions
% ===================================
u4 = y(:,1); % Vin (foot inlet)
p4 = y(:,2); % Pressure foot
u5 = y(:,3); % Uj, foot outlet, flue exit

% ===================================
% Interpolate computed solutions 
% (non-constant time step)
% at measured times. Errors calc.
% ===================================

F = griddedInterpolant(t,p4);
p4int = F(p3t);

relerror = abs(p4meas-p4int)./p4meas;
abserror = abs(p4meas-p4int);


% ===================================
%           Plot results
% ===================================

%%
figure(1); clf;

% /// velocity in ///
ax(1)=subplot(321);
plot(t,u4);ylabel('(m/s)');
title('(u4) Velocity foot-inlet [CALC]');

%  /// pressure foot ///
ax(2)=subplot(323);
plot(p3t,p4meas,'r');
hold on;
plot(t,p4,'b'); ylabel('(Pa)');
hold off;legend('measured','computed');
title('(p4) Pfoot (meas. vs simul.)');

% /// rel. error pressure foot ///
ax(3) = subplot(325);
plot(p3t,relerror,'o');ylabel('|meas-calc|/meas');
title('Pfoot: Relative error w.r.t. measured');
ylim([0. 1.2]);

% /// velocity out ///
ax(4)=subplot(322);
plot(t,u5); ylabel('(m/s)');
title('(u5) Jet Velocity [CALC]');
xlabel('time (s)');

% /// reduced jet velocity ///
ax(5)=subplot(324);
RDV = u5/(Wm*f1);
plot(t, RDV);
ylabel('\theta');
title('Reduced Jet Velocity [CALC]');

% /// abs. error pressure foot ///
ax(6) = subplot(326);
plot(p3t,abserror,'o');ylabel('(Pa)');
title('Pfoot: Abs. error Pfoot w.r.t. measured');



figure(2); clf;
%                      <<<<<  l5 EDIT >>>>>
ax(7) = plot(t(1:end-1), (2*(l5*0+5e-4)*diff(u5)./diff(t)) ./ u5(1:end-1).^2 ,'o');
ax(7) = ax(7).Parent;
title('2*l5*d(ujet)/dt  / ujet^2');
ylim(1*[-1,1]);

linkaxes(ax,'x');
xlim([0 5]);

function dydt = myodes(t,y,p3t,p3, pmmeas, params,mouthbool)

% Pass precomputed system coefficients
a1 = params(1);a2 = params(2);a3 = params(3);b1 = params(4);
b2 = params(5);c1 = params(6);c2 = params(7);c3 = params(8);

% Give measured P_resrv(t) pressure signal (p3) and 
% measuring time (p3t) to interpolate at integration time (t)
p3int = interp1(p3t,p3,t); 

dydt    = zeros(3,1);
dydt(1) = a1*p3int + a2*y(2) + a2*y(1)^2;
dydt(2) = b1*y(1) + b2*y(3);

if mouthbool % Account for measured mouth pressure
    p5int = interp1(p3t,pmmeas,t);
    
    % Amplitude oscillation in mouth 
    % w.r.t. pipe (Pp')->(Pm') (<Pm>=0)
        fac   = 1.1;  % 0.5 [TODO]: p_i Verge-Fourier
        p5int = fac*p5int;
        dydt(3) = c1*y(2) + c2*p5int + c3*y(3)^2;
    
    else % No mouth pressure (solve quickly)
    p5 = 0;
        dydt(3) = c1*y(2) + c2*p5 + c3*y(3)^2;
end
end

