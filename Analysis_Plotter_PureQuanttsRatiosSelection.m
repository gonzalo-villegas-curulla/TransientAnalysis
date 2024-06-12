
clc, clear;
try
    cd /media/organ/ExtremeSSD/OrganPipe2023-2024/DataTransients/processed
end

% ============= LOADS =========================

files = dir('*PROCESSED.mat');

% Resoantor length: (56)
load('./Geometry/Lp_m.mat'); 
LP = Lp_m; clear Lp_m

% Foot volume: (56)
load('./Geometry/Vf_m3.mat'); 
VF = Vf_m3; clear Vf_m3

% Width of pallet window slot (56)
load('./Geometry/KoenigPalletWindWidth_m.mat');
PW    = palletwinwidth_m; clear palletwinwidth_m

% Tone hole diam, on the wooden top board: (56)
load('./Geometry/ToneHoleDiam_m.mat');
TNHD = ToneHoleDiam; clear ToneHoleDiam

% Toe hole area (foot inlet): (56)
load('./Geometry/InletManipFitted_mm2.mat'); 
INLET = InletManipFitted_mm2* 1e-6; clear InletManipFitted_mm2

% Flue exit area (S_jet): (56)
load('./Geometry/FlueExitManipFitted_mm2.mat');  
SJET    = FlueExitManipFitted_mm2 * 1e-6; clear FlueExitManipFitted_mm2

% Flue exit little height: (56)
load('./Geometry/SmallhManipFitted_mm.mat');
h     = hManipFitted_mm * 1e-3; clear hManipFitted_mm

%FluExit width, big Hm: (56)
load('./Geometry/BigHManipFitted_mm.mat'); 
H     = HManipFitted_mm * 1e-3; clear HManipFitted_mm

% Mouth cutupt distance to labium (56):
load('./Geometry/WmManipFitted_mm.mat');
WM    = WmManipFitted_mm * 1e-3; clear WmManipFitted_mm 

% Resonator diameter (56)
load('./Geometry/DpManipFitted_mm.mat');
DP    = DpManipFitted_mm*1e-3; clear DpManipFitted_mm

% F1 (averaged) for the measured pipes (16) <<<<!
load('./Geometry/f1ManipMean.mat');
% >> F1MEAN

% ============= PARAMETERS =========================

rho   = 1.2;
c     = 340;
maskpipes  = [3,4,5,6,7,9,10,11,13,15,17,19,24,25,27, 29,32, 34,37,39,41,44]';



PalletDepth = 1e-3*129.8; %  [m]
PALLAREA    = PalletDepth * PW;  % Slot area covered by the valve [m^2]
VGROOVE     = PW *0.51*0.05;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefs and params 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vena contracti:
VCpallet = 0.62;
VCinlet  = 0.62;
VCjet    = 0.95;
VCpallet=1;VCinlet=1;VCjet=1;

NumTransMax = 61;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Q-factors for fundamentals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RADS     = 0.5*DP;
RADSmask = RADS(maskpipes);
DPmask   = DP(maskpipes);
LPmask   = LP(maskpipes);

n = 1; % mode
chi       = 3e-5 * sqrt(n*F1MEAN)./RADSmask;
OM        = n*pi*c./(LPmask + 1.2*RADSmask);
part1     = OM/(2*c);
part2     = (LPmask + 1.2*RADSmask)./(chi.*LPmask + (OM.^2.*RADSmask.^2)/(2*c^2) )  ;
QFAC1     = part1.*part2;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MX           = nan*ones(NumTransMax, length(files), 41);

vecmeanfreqs = zeros(length(files),1);
PRTMX        = nan*ones(NumTransMax, length(files), 3);
PRTpipeMedian = zeros(length(files),1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Populate data in matrices 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for idx = 1 : length(files)
    
    filename = files(idx).name;
    load(filename);
    
    numpipe    = str2num(filename(2:3));    
    
    % #56 data:
    oneLP      = LP(numpipe); % Pipe length
    oneVF      = VF(numpipe); % Volume foot 
    onePW      = PW(numpipe); % Pallet valve slot width 
    oneTNHD    = TNHD(numpipe);
    oneINLET   = INLET(numpipe);
    oneSJET    = SJET(numpipe);
    oneh       = h(numpipe);
    oneH       = H(numpipe); 
    oneWM      = WM(numpipe); % Mouth W cutup
    oneDP      = DP(numpipe); % Diameter pipe
    oneVGROOVE = VGROOVE(numpipe);
    
    % #16 data:
    onef1      = F1MEAN(idx);
    oneQFAC1   = QFAC1(idx);
    
    %%%%%%%%%% GEOMETRICAL ASPECTS %%%%%%%%%% 
%     NumTransMax = 1; % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    MX(:, idx, 1)    = oneLP*ones(NumTransMax,1);
    MX(:, idx, 2)    = oneVF*ones(NumTransMax,1);
    MX(:, idx, 3)    = onePW*ones(NumTransMax,1);
    MX(:, idx, 4)    = oneTNHD*ones(NumTransMax,1);
    MX(:, idx, 5)    = oneINLET*ones(NumTransMax,1);
    MX(:, idx, 6)    = oneSJET*ones(NumTransMax,1);
    MX(:, idx, 7)    = oneh*ones(NumTransMax,1);
    MX(:, idx, 8)    = oneH*ones(NumTransMax,1);
    MX(:, idx, 9)    = oneWM*ones(NumTransMax,1);
    MX(:, idx, 10)   = oneDP*ones(NumTransMax,1);
    MX(:, idx, 11)   = oneVGROOVE*ones(NumTransMax,1);
    MX(:, idx, 12)   = oneQFAC1*ones(NumTransMax,1);
    
    
    
    %%%%%%%%%% S-S ASPECTS %%%%%%%%%% 
    
    % S1/S2 factor, Spall/SJET(Pf,Pgrv) Remplissage factor:
    Remplissage = real(sqrt(Pfoot_targ ./ (Pgroove_targ-Pfoot_targ ) ) ); 
    Qpall2grv   = VCpallet * PALLAREA(numpipe) * real(sqrt( (PpalletB_targ - Pgroove_targ) *2/rho));
    Qgrv2ft     = VCinlet  * INLET(numpipe)    * real(sqrt( (Pgroove_targ  - Pfoot_targ  ) *2/rho));
    Qjet        = VCjet    * SJET(numpipe)     * real(sqrt( (Pfoot_targ    - 0           ) *2/rho));
    
    MX(find(f1)                 , idx, 13) = f1(:);
    MX(find(RJV)                , idx, 14) = RJV;
    MX(find(Qpall2grv),           idx, 15) = Qpall2grv;
    MX(find(Qgrv2ft)  ,           idx, 16) = Qgrv2ft;
    MX(find(Qjet),                idx, 17) = Qjet;
    MX(find(Remplissage)       ,  idx, 18) = Remplissage(Remplissage~=0);
    
    
    MX(find(Pfoot_targ), idx, 19) = PpalletB_targ;
    MX(find(Pfoot_targ), idx, 20) = Pgroove_targ;
    MX(find(Pfoot_targ), idx, 21) = Pfoot_targ;
    MX(find(Pfoot_targ), idx, 22) = Ppipe_targ;
    
    
    if 0
        MX(find(Pgroove_targ./PpalletB_targ), idx, 23)  = Pgroove_targ./PpalletB_targ;
        MX(find(Pfoot_targ./Pgroove_targ)   , idx, 24)  = Pfoot_targ  ./Pgroove_targ;
        MX(find(Ppipe_targ./Pfoot_targ)     , idx, 25)  = Ppipe_targ  ./Pfoot_targ;
    else
        MX(find(Pgroove_targ), idx, 23)     = Pgroove_targ - PpalletB_targ;
        MX(find(Pfoot_targ), idx, 24)       = Pfoot_targ   - Pgroove_targ;
        MX(find(Ppipe_targ), idx, 25)       = Ppipe_targ   - Pfoot_targ;
    end
    
    %%%%%%%%%% TRANSIENT ASPECTS %%%%%%%%%% 
    
    t20foot=t20foot(:);t20groove=t20groove(:);t20mouth=t20mouth(:);
    PRTgroove = PRTgroove(:);PRTfoot=PRTfoot(:);PRTpipe=PRTpipe(:);
    
    MX(find(betafit)  , idx, 26)  = betafit(:);
    MX(find(nufit)    , idx, 27)  = nufit(:);
    
    
    MX(find(PRTgroove), idx,28) = PRTgroove;
    MX(find(PRTfoot)  ,idx,29)  = PRTfoot;
    MX(find(PRTpipe)  ,idx,30)  = PRTpipe;
    
    MX(find(PRTfoot)     ,idx, 31) = PRTfoot./PRTgroove;
    MX(find(PRTpipe)     ,idx, 32) = PRTpipe./PRTfoot;
    
    MX(find(t20groove), idx, 33) = t20groove(t20groove~=0);
    MX(find(t20foot),   idx, 34) = t20foot(t20foot~=0);
    MX(find(t20mouth),  idx, 35) = t20mouth(t20mouth~=0);
    
    MX(find(t20foot-t20groove), idx, 36) = (t20foot-t20groove);
    MX(find(t20mouth-t20foot),  idx, 37) = (t20mouth-t20foot);
    
    MX(find(betafit)            , idx, 38)   = Area1(:);
    MX(find(betafit)            , idx, 39)   = Area2(:);
    
    
    MX(:           , idx, 40)   = ones(NumTransMax,1)*(INLET(numpipe)/PALLAREA(numpipe));
    MX(:           , idx, 41)   = ones(NumTransMax,1)*(SJET(numpipe)/INLET(numpipe));
    
    
   
    %%%%%%%%%% MARGINAL ASPECTS %%%%%%%%%% 
    PRTpipeMedian(idx) = median(PRTpipe(PRTpipe~=0),'omitnan');
    vecmeanfreqs(idx)  = mean(f1);
    freqref            = 440;
    freqlogax          = 12*log2(vecmeanfreqs/freqref);
    
    
end    
%%%%%%%%%%%%%%%%%%%% PREPARE PCA()  %%%%%%%%%%%%%%%%%%%%



MXresh = nan*ones(size(MX,1)*size(MX,2),size(MX,3));
for idx = 1 : size(MX,3)
   MXresh(:,idx) = reshape( MX(:,:,idx),numel(MX(:,:,idx)),1 ) ;  
end
    



% [1]:Lp        [2]:Vf       [3]:PWidth  [4]:ToneHoleDiam  [5]:Inlet 
% [6]:Sjet      [7]:h        [8]:H       [9]:Wm            [10]:Dpipe        
% [11]:Vgroove  [12]:Qfactor        

% [13]:f1          [14]:theta        [15]:Qpall2groove   
% [16]:Qgrv2foot   [17]:Qjet         [18]:Remplissage  
% [19]:Ppall targ  [20]:Pgrove targ  [21]:Pfoot targ    [22]:Prad targ

% [23]:Ptarg grv-pall     [24]:Ptarg foot-grv     [25]:Ptarg rad-foot

% [26]:beta             [27]:nu   
% [28]:PRTgrv           [29]:PRTfoot        [30]:PRTrad
% [31]:PRTfoot/grv      [32]:PRTrad/foot  
% [33]:t20 groove       [34]:t20 foot       [35]:t20 rad
% [36]:t20 foot-groove  [37]:t20 rad-foot
% [38]:Area1            [39]:Area2
% [40]:Sin/Spall [41]:Sjet/Sin


% maskpca = [16,17,23,24,25,28,29,30,33,34,35,38,39];
maskpca = [16,17,19:22,28,29,30,33,34,35,38,39];
% maskpca = [16,17,24,29,30,33,34,35,38,39];

namevarsall = {'Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'P{o} pall','P{o}grv','P{o}foot','P{o}rad',...
    'P{o} grv-pall','P{o}foot-grv','P{o}rad-foot',...
    'Beta','nu',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't{20}grv','t{20}foot','t{20}rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin'};

namevars2 = ['Lp','Vf','PWdth','TNHD','INLET','SJET','h','H','Wm','Dp',...
    'Vgrv','Qfact','f1','theta','Qpall2grv','Qgr2ft','Qjet','Rmplssg',...
    'PpallTarg','PgrvTarg','PfootTarg','PradTarg',...
    'Pgrv-Ppall(TARG)','Pfoot-Pgrv(TARG)','Prad-Pfoot(TARG)',...
    'PRTgrv','PRTfoot','PRTrad',...
    'PRTfoot/grv','PRTrad/foot',...
    't20grv','t20foot','t20rad',...
    't20(foot-grv)','t20(rad-foot)',...
    'Area1','Area2','Sin/Spal','Sjet/Sin'];

namevars3 = ["Lp","Vf","PWdth","TNHD","INLET","SJET","h","H","Wm","Dp",...
    "Vgrv","Qfact","f1","theta","Qpall2grv","Qgr2ft","Qjet","Rmplssg",...
    "PpallTarg","PgrvTarg","PfootTarg","PradTarg",...
    "Pgrv-Ppall(TARG)","Pfoot-Pgrv(TARG)","Prad-Pfoot(TARG)",...
    "PRTgrv","PRTfoot","PRTrad",...
    "PRTfoot/grv","PRTrad/foot",...
    "t20grv","t20foot","t20rad",...
    "t20(foot-grv)","t20(rad-foot)",...
    "Area1","Area2","Sin/Spal","Sjet/Sin"];

namevars4 = {"Lp","Vf","PWdth","TNHD","INLET","SJET","h","H","Wm","Dp",...
    "Vgrv","Qfact","f1","theta","Qpall2grv","Qgr2ft","Qjet","Rmplssg",...
    "P^{o}pall","P^{o}grv","P^{o}foot","P^{o}rad",...
    "P^{o}grv-Ppall","P^{o}foot-Pgrv","P^{o}rad-foot",...
    "PRTgrv","PRTfoot","PRTrad",...
    "PRTfoot/grv","PRTrad/foot",...
    "t20grv","t20foot","t20rad",...
    "t20(foot-grv)","t20(rad-foot)",...
    "Area1","Area2","Sin/Spal","Sjet/Sin"};



% With respect to beta
idxanalys = 26;
selmask   = [1:14]; % Only geometrical parameters 
namesmask = namevars4(selmask);
if 1

%     [idx, scores] = fsrftest(MXresh(:,[1:10,12:25,28:39]), MXresh(:,idxanalys));
    
    [idx, scores] = fsrftest(MXresh(:,selmask), MXresh(:,idxanalys));
    figure(1);clf;
    bar(scores(idx));
    set(gca, 'xtick',[1:length(idx)]);
    set(gca, 'xticklabels', namesmask(idx));
    title(sprintf('For: %s',namevars4{idxanalys}));
% else
    mdl = fsrnca(MXresh(:,selmask), MXresh(:,idxanalys));
    figure(2);clf;
    plot(mdl.FeatureWeights,'ro');
    grid on;
    xlabel('Feature index');
    ylabel('Feature weight');
    IDX = find( mdl.FeatureWeights > (mean(mdl.FeatureWeights)) );
    hold on;
    for jdx =1:length(IDX)
        text(IDX(jdx),mdl.FeatureWeights(IDX(jdx)),namesmask{IDX(jdx)});
    end
    xlim([0,length(mdl.FeatureWeights)+1]);
    
end 
% % % [idx, scores] = fsrmrmr(MXresh(:,[1:10,12:25,28:39]), MXresh(:,26));




%% PRT*f1
figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,13).*MX(:,:,29) ), 'b', 'filled');
ylabel('PRT * f_1 (ln)');
xlabel('$12log_2(f_1/440)$');

%% PRT rad / PRT foot
figure();
scatter( 12*log2(MX(:,:,13)/440), log(MX(:,:,32) ), 'b', 'filled');
ylabel('PRTrad / PRT_{foot} (ln)');
xlabel('$12log_2(f_1/440)$');
%% PRT foot / PRT groove
figure();
scatter( 12*log2(MX(:,:,13)/440), log(MX(:,:,31) ), 'b', 'filled');
ylabel('PRT$_{foot}$ / PRT$_{grv}$ (ln)');
xlabel('$12log_2(f_1/440)$');

%% t20foot / t20 groove

figure();
scatter( 12*log2(MX(:,:,13)/440),log(MX(:,:,34)./MX(:,:,33) ), 'b', 'filled');
ylabel('t$^{20}_{foot}$ / t$^{20}_{grv}$ (ln)');
xlabel('$12log_2(f_1/440)$');

%% t20rad / t20 foot

figure();
scatter( 12*log2(MX(:,:,13)/440),log(MX(:,:,35)-MX(:,:,34) ), 'b', 'filled');
ylabel('t$^{20}_{rad}$ - t$^{20}_{foot}$ (ln)');
xlabel('$12log_2(f_1/440)$');

%% Sj / Sin 

figure();
scatter( 12*log2(MX(:,:,13)/440), log(MX(:,:,6)./MX(:,:,5) ), 'b', 'filled');
ylabel('$\mathcal{S}_j / \mathcal{S}_{in}$ (ln)');
xlabel('$12log_2(f_1/440)$');

%% Pgrv / Ppall

figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,20)./MX(:,:,19) ), 'b', 'filled');
ylabel('$P_{grv}/P_{pall}$');
xlabel('$12log_2(f_1/440)$');

%% Pfoot / Pgrv

figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,21)./MX(:,:,20) ), 'b', 'filled');
ylabel('$P_{foot}/P_{grv}$');
xlabel('$12log_2(f_1/440)$');

%% Prad / Pfoot

figure();
scatter( 12*log2(MX(:,:,13)/440), (MX(:,:,22)./MX(:,:,21) ), 'b', 'filled');
ylabel('$P^o_{rad}/P^o_{foot}$');
xlabel('$12log_2(f_1/440)$');

%% t80 foot

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,29)+MX(:,:,34) ), 'b', 'filled');
ylabel('$t^{80}_{foot}$ [ms]');
xlabel('$12log_2(f_1/440)$');

%% t80 foot / T1 = t80foot *f1

figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,13).*(MX(:,:,29)+MX(:,:,34) ), 'b', 'filled');
ylabel('$t^{80}_{foot} * f_1$');
xlabel('$12log_2(f_1/440)$');


%% t80 rad

figure();
scatter( 12*log2(MX(:,:,13)/440), 1e3*(MX(:,:,30)+MX(:,:,35) ), 'b', 'filled');
ylabel('$t^{80}_{rad}$ [ms]');
xlabel('$12log_2(f_1/440)$');

%% t80 rad / T1 = t80rad *f1

figure();
scatter( 12*log2(MX(:,:,13)/440), MX(:,:,13).*(MX(:,:,30)+MX(:,:,35) ), 'b', 'filled');
ylabel('$t^{80}_{rad} * f_1$');
xlabel('$12log_2(f_1/440)$');




%% VC factor groove

factor = (MX(:,:,5)./(MX(:,:,3)*PalletDepth)).*sqrt( (MX(:,:,20) - MX(:,:,21) )./(MX(:,:,19)-MX(:,:,20)));

figure();
scatter( MX(:,:,5)*1e6, factor, 'b', 'filled');box on;
xlabel('$S^{geo}_{in} \ [mm^2]$','interpreter','latex','fontsize',FSZ);
ylabel('$\Gamma_{grv}$','interpreter','latex','fontsize',FSZ);





%% #################################################################### 
% Ratio of effective areas as a function of geometric areas ratio
% #################################################################### 


FSZ = 22;
figure();

% S foot/pallet ====

datax = (MX(:,:,5))./(MX(:,:,3)*PalletDepth); % Sin/Sgroove GEOMETRIC section
datay = sqrt( (  MX(:,:,19)-MX(:,:,20)  )./(MX(:,:,20)-MX(:,:,21)) ); % Effective section (SS)

axeshandle(1) = subplot(121);
scatter( datax , datay , 'b', 'filled');box on;
xlabel('$\frac{S^{\ geo}_{foot}}{S^{\ geo}_{pall}}$','interpreter','latex','fontsize',FSZ+8);
ylabel('$\frac{S^{\ eff}_{foot}}{S^{\ eff}_{pall}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);
set( get(gca,'YAxis'), 'FontSize', FSZ-8);

% S jet/foot ====

dataone = MX(:,:,6)./MX(:,:,5); % GEOMETRIC section Sjet/Sin,foot
datatwo = sqrt( (MX(:,:,20)-MX(:,:,21)   )./( MX(:,:,21)  ));


axeshandle(2) = subplot(122);
scatter( dataone, datatwo, 'b', 'filled');box on;
hold on;
plot([0,1.4],[0,1.4],'--r');
xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);
set( get(gca,'YAxis'), 'FontSize', FSZ-8);

%%
% #################################################################### 
% GOOD physical meaning and expected values
% #################################################################### 

%(VC factor foot)

FSZ = 22;
dataone = MX(:,:,6)./MX(:,:,5); % GEOMETRIC section Sjet/Sin,foot
datatwo = sqrt( (MX(:,:,20)-MX(:,:,21)   )./( MX(:,:,21)  ));

figure();
scatter( MX(:,:,5)*1e6, dataone./datatwo, 'b', 'filled');box on;
% scatter( dataone, dataone./datatwo, 'b', 'filled');box on;
% scatter( sqrt(MX(:,:,5)/pi)./MX(:,:,3), dataone./datatwo, 'b', 'filled');box on;
xlabel('$S^{geo}_{in} \ [mm^2]$','interpreter','latex','fontsize',FSZ);
ylabel('$factor$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
% grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ);
set( get(gca,'YAxis'), 'FontSize', FSZ);
ylim([0,1]);
%%
figure(); 
% scatter(MX(:,:,3), MX(:,:,5), 'b', 'filled');
plot( (MX(:,:,3)*0.05./MX(:,:,5))' , 'o');

%%

figure();
scatter( MX(:,:,5), datatwo, 'b', 'filled');box on;
xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ);

%%
XX = datax .* dataone;
YY = datay .* datatwo;

figure();

scatter( XX, YY, 'b', 'filled');box on;
% xlabel('$\frac{S^{geo}_{jet}}{S^{geo}_{foot}}$','interpreter','latex','fontsize',FSZ);
% ylabel('$\frac{S^{eff}_{jet}}{S^{eff}_{foot}}$','interpreter','latex','fontsize',FSZ);
set(get(gca,'ylabel'),'Rotation',0);
grid on; ylim([0 2]);
set( get(gca,'XAxis'), 'FontSize', FSZ-8);


%% 
% #################################################################### 
% AUTO CORRELATION OF VARIABLE VALUES 
% #################################################################### 
tmp = MXresh;

for idx = 1 : size(tmp,1)
    for jdx = 1 : size(tmp,2)
        if isnan(tmp(idx,jdx))
            tmp(idx,jdx) = 0;
        end
    end
end

DataToWorkWith = MXresh(:,maskpca);


tmp = (tmp-mean(tmp,1,'omitnan'))./var(tmp,'omitnan');

% tmp = log(abs(tmp)+1e-6);
tmp = exp(tmp);



tmpc = corrcoef(tmp);

tmpc(abs(tmpc)<0.0001) = 0;

figure();
imagesc(tmpc);
ax=gca;
ax.YDir = 'normal';
ax.XTick = 1:length(tmpc);
ax.XTickLabel = namevarsall;
ax.YTick = 1:length(tmpc);
ax.YTickLabel = namevarsall;
% colormap(bluewhitered); 

colorbar;


%%
% #################################################################### 
% PCA()
% #################################################################### 

% maskpca = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39];
% maskpca = [2,3,5,6,7,26,38];


DataToWorkWith = MXresh(:,maskpca);


wei = (DataToWorkWith-mean(DataToWorkWith,1,'omitnan'))./var(DataToWorkWith,'omitnan');

wei = ones(length(maskpca),1)';

[coeff,score,latent,tsquared,explained] = pca(DataToWorkWith, 'VariableWeights',wei);
figure(3);clf; plot(explained,'o-');
LBL = append(namevarsall(maskpca));


figure(4);clf;
biplot(coeff(:,1:3),'scores',score(:,1:3),'VarLabels',LBL);

%% 
% #################################################################### 
% BETA and NU with FREQUENCY
% #################################################################### 

figure(25); clf;
FSZ = 15;
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;

ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);

%% 
% #################################################################### 
% BETA and NU with FREQUENCY
% #################################################################### 

figure(25); clf;
FSZ = 15;
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( freqlogax, log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;

ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);

%% 

% #################################################################### 
% BETA and NU with I21 and I31
% #################################################################### 
figure(); clf;
FSZ = 15;
scatter( log(MX(:,:,38)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
xlabel('I21');

figure();
scatter( log(MX(:,:,38)), log(MX(:,:,27)), 'r', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('I21');

%% 
% #################################################################### 
% BETA and NU with PTARG GROOVE
% #################################################################### 


figure(27); clf;
FSZ = 15;
scatter( log(MX(:,:,20)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log(MX(:,:,20)), log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log PTARG groove$','interpreter','latex','fontsize',FSZ);



%% 
% #################################################################### 
% BETA and NU with PTARG FOOT
% #################################################################### 


figure(28); clf;
FSZ = 15;
scatter( log(MX(:,:,21)), log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
hold on;
yyaxis right;
scatter( log(MX(:,:,21)), log(MX(:,:,27)), 'r', 'filled');box on;
% ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
% legend('$\beta$','$\nu$','interpreter','latex','fontsize',FSZ);
yyaxis right;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
xlabel('$log PTARG foot$','interpreter','latex','fontsize',FSZ);


%% 

% #################################################################### 
% PRESSURE DROPS
% #################################################################### 

FSZ = 14;

fili = [1:9,12,15:22];
fili = [1:22];
DROP1 = MX(:,fili,20)-MX(:,fili,19) ;
DROP2 = MX(:,fili,21)-MX(:,fili,20);
DROP3 = -MX(:,fili,21);%-MX(:,fili,19);

if 0
figure(24);clf;
sbh(1) = subplot(1,2,1);
sasa=pcolor(DROP1); colormap(inferno);colorbar; sasa.FaceColor = 'interp'; sasa.LineStyle = 'none';
title(sprintf('Pallet2Groove. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP1,[1,2],'omitnan'), std(DROP1,0,[1 2],'omitnan')   ));

sbh(2) = subplot(1,2,2);
soso=pcolor(DROP2); colormap(inferno);colorbar; soso.FaceColor = 'interp'; soso.LineStyle = 'none';
title(sprintf('Groove2Foot. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP2,[1,2],'omitnan'), std(DROP2,0,[1 2],'omitnan')   ));

linkaxes(sbh,'xy');
ylim([1 36])

figure(13); clf; 
sisi = pcolor(DROP3);colormap(inferno);colorbar;sisi.FaceColor='interp';sisi.LineStyle='none';ylim([1 38]);
title(sprintf('Pallet2Foot. Mean drop: %1.2f [Pa], std %1.3f [Pa]',mean(DROP3,[1,2],'omitnan'), std(DROP3,0,[1 2],'omitnan')   ));
end

figure(28);clf;
scatter(fili, DROP1,'b');
hold on;
scatter(fili, DROP2, 'r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);
ylabel('Pressure drop [Pa]','interpreter','latex','fontsize',FSZ);
scatter(fili, DROP3, 'g');







%% 
% #################################################################### 
% VENA CONTRACTA Inlet-Jet (steady state values)
% #################################################################### 

FSZ = 16;



MXh   = MX(:,:,7);
MXH   = MX(:,:,8);
MXpal = MX(:,:,19);
MXgr  = MX(:,:,20);
MXft  = MX(:,:,21);

rat2 = (MXh.*MXH)./(MX(:,:,5)).*sqrt( (MXft)./(MXgr-MXft)  );

figure();
RatioToPlot = 1./rat2;
linidx = 1:min(size(RatioToPlot));
scatter(linidx, RatioToPlot, 'b','filled');%ylim([0 1]);
title('Vena contracta ratios: $VC_{jet}/VC_{inlet}$','interpreter','latex','fontsize',FSZ);
p = polyfit(linidx, mean(RatioToPlot,1,'omitnan'),1);
hold on;
plot(polyval(p,linidx),'--r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);


%% 

% #################################################################### 
% VENA CONTRACTA palletbox-groove-foot
% #################################################################### 


PALLAREA = repmat(PALLAREA(maskpipes)',max(size(MXgr)),1);

rat1 = MX(:,:,5)./PALLAREA .*sqrt( (MXgr-MXft) ./ (MXpal-MXgr) );

figure();
RatioToPlot = 1./rat1;
linidx = 1:min(size(RatioToPlot));
scatter(linidx, RatioToPlot, 'b','filled');%ylim([0 1]);
title('Vena contracta ratios: $VC_{inlet}/VC_{groove-slot}$','interpreter','latex','fontsize',FSZ);
p = polyfit(linidx, mean(RatioToPlot,1,'omitnan'),1);
hold on;
plot(polyval(p,linidx),'--r');
xlabel('Num pipe','interpreter','latex','fontsize',FSZ);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 17;
figure();


geoh(1) = subplot(3,4,1);
scatter( freqlogax, log(MX(:,:,1)), 'b', 'filled');box on;
ylabel('$Lp$','interpreter','latex','fontsize',FSZ);

geoh(2) = subplot(3,4,2);
scatter( freqlogax, log(MX(:,:,2)), 'b', 'filled');box on;
ylabel('$Vf$','interpreter','latex','fontsize',FSZ);

geoh(3) = subplot(3,4,3);
scatter( freqlogax, log(MX(:,:,3)), 'b', 'filled');box on;
ylabel('$PWdth$','interpreter','latex','fontsize',FSZ);

geoh(4) = subplot(3,4,4);
scatter( freqlogax, log(MX(:,:,4)), 'b', 'filled');box on;
ylabel('$TnHD$','interpreter','latex','fontsize',FSZ);

geoh(5) = subplot(3,4,5);
scatter( freqlogax, log(MX(:,:,5)), 'b', 'filled');box on;
ylabel('$Sin$','interpreter','latex','fontsize',FSZ);

geoh(6) = subplot(3,4,6);
scatter( freqlogax, log(MX(:,:,6)), 'b', 'filled');box on;
ylabel('$Sjet$','interpreter','latex','fontsize',FSZ);

geoh(7) = subplot(3,4,7);
scatter( freqlogax, log(MX(:,:,7)), 'b', 'filled');box on;
ylabel('$h$','interpreter','latex','fontsize',FSZ);

geoh(8) = subplot(3,4,8);
scatter( freqlogax, log(MX(:,:,8)), 'b', 'filled');box on;
ylabel('$H$','interpreter','latex','fontsize',FSZ);

geoh(9) = subplot(3,4,9);
scatter( freqlogax, log(MX(:,:,9)), 'b', 'filled');box on;
ylabel('$Wm$','interpreter','latex','fontsize',FSZ);

geoh(10) = subplot(3,4,10);
scatter( freqlogax, log(MX(:,:,10)), 'b', 'filled');box on;
ylabel('$Dp$','interpreter','latex','fontsize',FSZ);

% geoh(11) = subplot(3,4,11);
% scatter( freqlogax,log( MX(:,:,11)), 'b', 'filled');box on;
% ylabel('$Vgrv$','interpreter','latex','fontsize',FSZ);

geoh(12) = subplot(3,4,12);
scatter( freqlogax, log(MX(:,:,12)), 'b', 'filled');box on;
ylabel('$Qfact1$','interpreter','latex','fontsize',FSZ);


xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);
linkaxes(geoh,'x');    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 2-3-4 (S-S) (pall-groove-foot)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 17;

figure();
f1h(1)=subplot(8,2,1);
scatter( freqlogax, MX(:,:,15), 'b', 'filled');box on;
ylabel('$Q_{pall2gr}$','interpreter','latex','fontsize',FSZ); title('Steady-State');

f1h(2)=subplot(8,2,2); 
scatter( freqlogax, log(MX(:,:,26)), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient');
f1h(3)=subplot(8,2,3);
scatter( freqlogax, MX(:,:,16), 'b', 'filled');box on;
ylabel('$Q_{gr2ft}$','interpreter','latex','fontsize',FSZ);
f1h(4)=subplot(8,2,4);
scatter( freqlogax, MX(:,:,27), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);
f1h(5)=subplot(8,2,5);
scatter( freqlogax, MX(:,:,19), 'b', 'filled');box on;
ylabel('$P^{o}_{pall}$','interpreter','latex','fontsize',FSZ);
f1h(6)=subplot(8,2,8);
scatter( freqlogax, MX(:,:,28), 'b', 'filled');box on;
ylabel('$PRT_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(7)=subplot(8,2,7);
scatter( freqlogax, MX(:,:,20), 'b', 'filled');box on;
ylabel('$P^{o}_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(8)=subplot(8,2,10);
scatter( freqlogax, MX(:,:,29), 'b', 'filled');box on;
ylabel('$PRT_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(9)=subplot(8,2,9);
scatter( freqlogax, MX(:,:,21), 'b', 'filled');box on;
ylabel('$P^{o}_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(10)=subplot(8,2,12);
scatter( freqlogax, MX(:,:,31), 'b', 'filled');box on;
ylabel('$PRT_{(foot/grv)}$','interpreter','latex','fontsize',FSZ);
f1h(11)=subplot(8,2,11);
scatter( freqlogax, MX(:,:,23), 'b', 'filled');box on;
ylabel('$P^{o}_{(grv-pall)}$','interpreter','latex','fontsize',FSZ);
f1h(12)=subplot(8,2,14);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;ylim([0 0.010]);
ylabel('$t^{20}_{grv}$','interpreter','latex','fontsize',FSZ);
f1h(13)=subplot(8,2,13);
scatter( freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$P^{o}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);
% f1h(14)=subplot(8,2,14);
% scatter( freqlogax, MX(:,:,34), 'b', 'filled');box on;ylim([0 0.010]);
% ylabel('$t^{20}_{foot}$','interpreter','latex','fontsize',FSZ);
f1h(15)=subplot(8,2,16);
scatter( freqlogax, MX(:,:,36), 'b', 'filled');box on;ylim([0 0.003]);
ylabel('$t^{20}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 3-4-5 (S-S) (groove-foot-outside)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FSZ = 14;

figure();
f2h(1)=subplot(8,3,19);
scatter( freqlogax, MX(:,:,13), 'b', 'filled');box on;
ylabel('$f_1$','interpreter','latex','fontsize',FSZ);

f2h(2)=subplot(8,3,2);
scatter( freqlogax, MX(:,:,26), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);title('Transient')

f2h(3)=subplot(8,3,6);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;
ylabel('$t^{20}_{grv}$','interpreter','latex','fontsize',FSZ);ylim([0 10e-3]);

f2h(4)=subplot(8,3,1);
scatter( freqlogax, MX(:,:,14), 'b', 'filled');box on;ylim([0 10]);
ylabel('$\theta$','interpreter','latex','fontsize',FSZ);title('Steady-State');

f2h(5)=subplot(8,3,3);
scatter( freqlogax, MX(:,:,27), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ);title('Transient')

f2h(6)=subplot(8,3,9);
scatter( freqlogax, MX(:,:,34), 'b', 'filled');box on;
ylabel('$t^{20}_{foot}$','interpreter','latex','fontsize',FSZ);ylim([0 10e-3]);

f2h(7)=subplot(8,3,22);
scatter( freqlogax, MX(:,:,17), 'b', 'filled');box on;
ylabel('$Q_{jet}$','interpreter','latex','fontsize',FSZ);

f2h(8)=subplot(8,3,5);
scatter( freqlogax, MX(:,:,28), 'b', 'filled');box on;
ylabel('$PRT_{grv}$','interpreter','latex','fontsize',FSZ);

f2h(9)=subplot(8,3,12);
scatter( freqlogax, MX(:,:,35), 'b', 'filled');box on;ylim([0 0.03])
ylabel('$t^{20}_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(10)=subplot(8,3,4);
scatter( freqlogax, MX(:,:,20), 'b', 'filled');box on;
ylabel('$P^{o}_{grv}$','interpreter','latex','fontsize',FSZ);

f2h(11)=subplot(8,3,8);
scatter( freqlogax, MX(:,:,29), 'b', 'filled');box on;
ylabel('$PRT_{foot}$','interpreter','latex','fontsize',FSZ);

f2h(12)=subplot(8,3,15);
scatter( freqlogax, MX(:,:,36), 'b', 'filled');box on;ylim([0 0.003]);
ylabel('$t^{20}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);

f2h(13)=subplot(8,3,7);
scatter( freqlogax, MX(:,:,21), 'b', 'filled');box on;
ylabel('$P^{o}_{foot}$','interpreter','latex','fontsize',FSZ);

f2h(14)=subplot(8,3,11);
scatter( freqlogax, MX(:,:,30), 'b', 'filled');box on;
ylabel('$PRT_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(15)=subplot(8,3,18);
scatter( freqlogax, MX(:,:,37), 'b', 'filled');box on;ylim([0 0.02]);
ylabel('$t^{20}_{(rad-ft)}$','interpreter','latex','fontsize',FSZ);

f2h(16)=subplot(8,3,10);
scatter( freqlogax, MX(:,:,22), 'b', 'filled');box on;
ylabel('$P^{o}_{rad}$','interpreter','latex','fontsize',FSZ);

f2h(17)=subplot(8,3,14);
scatter( freqlogax, MX(:,:,31), 'b', 'filled');box on;
ylabel('$PRT_{(ft/grv)}$','interpreter','latex','fontsize',FSZ);

f2h(18)=subplot(8,3,[20,23]);
scatter( freqlogax, log(MX(:,:,38)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I1 \ (log)$','interpreter','latex','fontsize',FSZ);

f2h(19)=subplot(8,3,13);
scatter( freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$P^{o}_{(ft-grv)}$','interpreter','latex','fontsize',FSZ);

f2h(20)=subplot(8,3,17);
scatter( freqlogax, MX(:,:,32), 'b', 'filled');box on;
ylabel('$PRT_{(rad/ft)}$','interpreter','latex','fontsize',FSZ);

f2h(21)=subplot(8,3,[21,24]);
scatter( freqlogax, log(MX(:,:,39)), 'b', 'filled');box on;ylim([-9 -1]);
ylabel('$I2 \ (log)$','interpreter','latex','fontsize',FSZ);

f2h(22)=subplot(8,3,16);
scatter( freqlogax, MX(:,:,25), 'b', 'filled');box on;
ylabel('$P^{o}_{(rad-ft)}$','interpreter','latex','fontsize',FSZ);


linkaxes(f2h,'x');
xlim([-21 23]);

    
%% %%%% PLOT %%%%%
    
% BENOIT'S PLOT (not the correct layers any more)
FSZ = 16;

figure(5); clf;
Bh(1) = nexttile([1 12]);
scatter( freqlogax, MX(:,:,33), 'b', 'filled');box on;
ylabel('$P^o_{foot}$','interpreter','latex','fontsize',FSZ);

Bh(2) = nexttile([1 12]);
scatter(freqlogax,  MX(:,:,6)./ MX(:,:,5), 'b', 'filled');box on;
ylabel('$S_{jet}/S_{in}$','interpreter','latex','fontsize',FSZ);

Bh(3) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,16), 'b', 'filled');box on;
ylabel('$Q_{in}$','interpreter','latex','fontsize',FSZ);

Bh(4) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,14), 'b', 'filled');box on;
ylabel('$\theta$','interpreter','latex','fontsize',FSZ);

Bh(5) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,24), 'b', 'filled');box on;
ylabel('$GRV_{PRT}$','interpreter','latex','fontsize',FSZ);

Bh(6) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,25), 'b', 'filled');box on;
ylabel('$FOOT_{PRT}$','interpreter','latex','fontsize',FSZ);

Bh(7) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,23), 'b', 'filled');box on;
ylabel('$\nu$','interpreter','latex','fontsize',FSZ)

Bh(8) = nexttile([1 12]);
scatter(freqlogax, MX(:,:,22), 'b', 'filled');box on;
ylabel('$\beta$','interpreter','latex','fontsize',FSZ);

xlabel('$12log_2(f_1 / f_{ref})$','interpreter','latex','fontsize',FSZ);
linkaxes(Bh,'x');    
    
    
    
    
    
    
    
    
    
    
    
    
%%
    


scatter(axh(1),freqlogax, MX(:,:,1),'filled','b');
scatter(axh(2),freqlogax, MX(:,:,2),'filled','b');
scatter(axh(3),freqlogax, MX(:,:,3),'filled','b');

scatter(axh(4),freqlogax, MX(:,:,4),'filled','b');

scatter(axh(5),freqlogax, MX(:,:,5),'filled','b');
scatter(axh(6),freqlogax, MX(:,:,6),'filled','b');

scatter(axh(7),freqlogax, MX(:,:,7),'filled','b');ylim(axh(7),[0.5 1]);
scatter(axh(8),freqlogax, MX(:,:,8),'filled','b');ylim(axh(8),[0 1]);
scatter(axh(9),freqlogax, MX(:,:,9),'filled','b'); ylim(axh(9), [0 0.25]);
scatter(axh(10),freqlogax, MX(:,:,10),'filled','b');%ylim(axh(10),[0 3]);

scatter(axh(11),freqlogax, MX(:,:,11),'filled','b');%ylim(axh(11),[-0 0.4]);
scatter(axh(12),freqlogax, MX(:,:,12),'filled','b');ylim(axh(12),[1.3 2.1]);
scatter(axh(13),freqlogax, MX(:,:,13),'filled','b');ylim(axh(13),[0 30]);
scatter(axh(14),freqlogax, MX(:,:,14),'filled','b');
scatter(axh(15),freqlogax, MX(:,:,15),'filled','b');

scatter(axh(16),freqlogax,  MX(:,:,16),'filled','b');
scatter(axh(17),freqlogax,  MX(:,:,17),'filled','b');
scatter(axh(18),freqlogax,  MX(:,:,18),'filled','b');

scatter(axh(19),freqlogax, MX(:,:,19),'filled','b');

linkaxes(axh(:),'x');
xlim([0.99*min(freqlogax),1.05*max(freqlogax)]);
grid(axh, 'on');box(axh, 'on');
drawnow();


xlim([-20 25]);


%% 

% #################################################################### 
% PRT alones 
% #################################################################### 

figure();
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,1),'filled','b');%ylim([0.5 1]);
ylabel('groove'); ylim([0 0.005]);
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,2),'filled','b');%ylim([0.5 1]);
ylabel('foot');ylim([0 0.005]);
nexttile([1,4]);   
scatter(freqlogax, PRTMX(:,:,3),'filled','b');%ylim([0.5 1]);
ylabel('pipe');ylim([0 0.200]);
hold on;
plot(freqlogax, PRTpipeMedian,'ro-');


%% 

% #################################################################### 
% Q-factor vs BETA
% #################################################################### 

figure(12);clf; hold on;
% Q-factor vs Beta
try
for idx = 1 : length(QFAC1)
    scatter(QFAC1(idx), log(MX(:, idx, 5)), 'b','filled');
end
end
xlabel('$Q-$factor [s$^{-1}$]','interpreter','latex','fontsize',14);
ylabel('$\beta$ fitted (log) [s$^{-1}$]','interpreter','latex','fontsize',14);
hold on

BM = MX(:, :, 5);
BM = mean(BM,1, 'omitnan');
p = polyfit(QFAC1(find(BM)), log(BM), 1);
betalin = polyval(p, QFAC1(find(BM)) );
plot(QFAC1(find(betalin)), (betalin), '--r');

box on;
%% 

% #################################################################### 
% NU versus Sj/Sin
% #################################################################### 

figure(13); clf; hold on;

palletww_mask   = PW(maskpipes);
palletarea_mask = PALLAREA(maskpipes);
Vf_mask         = VF(maskpipes);
hmask           = h(maskpipes);
Wmmask          = WM(maskpipes);
Inletmask       = INLET(maskpipes);
Sjetmask        = SJET(maskpipes);

NU = MX(:, :, 27);
ratios = Sjetmask./Inletmask;
try
for idx = 1 : size(NU,2)
scatter(  log(Vf_mask(idx)), log(NU(:,idx)), 'b','filled'); % YES good
% scatter( log( ratios(idx) ), log((NU(:,idx))), 'b','filled'); 
% scatter( log( Inletmask(idx)./palletarea_mask(idx) ), log((NU(:,idx))),'b','filled'); 
% scatter( log( Sjetmask(idx)./palletarea_mask(idx) ), log((NU(:,idx))),'b','filled'); 
% scatter( log( Inletmask(idx) ), log((NU(:,idx))),'b','filled'); 

end
end


hold on;

NM = mean(NU,1, 'omitnan');
p = polyfit(log(ratios), log(NM), 1);
nulin = polyval(p, log(ratios) );
% plot(log(ratios), (nulin), '--r')


FSZ = 14;
xlabel('$S_{jet}/S_{inlet}$ (log)','interpreter','latex','fontsize',FSZ);
ylabel('$\nu$ (ln)','interpreter','latex','fontsize',FSZ);
box on;