function [Lp,Vf,Pw,Tnhd,Sin,Sj,Hm,h,Wm,Rp,palletLHS,palletWid,palletRHS,palletHtraj] = getgeometry(PIPENUM)
pipenames = [03,04,05,06,07,09,10,11,13,15,17,19,24,25,27,29,32,34,38,39,41,44 ];


load('../Geometry/Lp_m.mat'); 
load('../Geometry/Vf_m3.mat'); 
load('../Geometry/KoenigPalletWindWidth_m.mat');
load('../Geometry/ToneHoleDiam_m.mat');
load('../Geometry/InletManipFitted_mm2.mat'); 
load('../Geometry/FlueExitManipFitted_mm2.mat');  
load('../Geometry/SmallhManipFitted_mm.mat');
load('../Geometry/BigHManipFitted_mm.mat'); 
load('../Geometry/WmManipFitted_mm.mat');
load('../Geometry/DpManipFitted_mm.mat');

LP    = Lp_m; clear Lp_m
VF    = Vf_m3; clear Vf_m3
PW    = palletwinwidth_m; clear palletwinwidth_m
TNHD  = ToneHoleDiam; clear ToneHoleDiam
SIN   = InletManipFitted_mm2* 1e-6; clear InletManipFitted_mm2
SJET  = FlueExitManipFitted_mm2 * 1e-6; clear FlueExitManipFitted_mm2
HM    = HManipFitted_mm * 1e-3; clear HManipFitted_mm
WM    = WmManipFitted_mm * 1e-3; clear WmManipFitted_mm 
DP    = DpManipFitted_mm*1e-3; clear DpManipFitted_mm
h     = hManipFitted_mm * 1e-3; clear hManipFitted_mm


Lp   = LP(pipenames(PIPENUM));
Vf   = VF(pipenames(PIPENUM));
Pw   = PW(pipenames(PIPENUM));
Tnhd = TNHD(pipenames(PIPENUM));
Sin  = SIN(pipenames(PIPENUM));
Sj   = SJET(pipenames(PIPENUM));
Hm   = HM(pipenames(PIPENUM));
h    = h(pipenames(PIPENUM));
Wm   = WM(pipenames(PIPENUM));
Dp   = DP(pipenames(PIPENUM));
Rp   = 0.5*Dp;

run PalletValveDimensions.m

palletLHS   = PalletValveGeometry(pipenames(PIPENUM),1);
palletWid   = PalletValveGeometry(pipenames(PIPENUM),2);
palletRHS   = PalletValveGeometry(pipenames(PIPENUM),3);
palletHtraj = PalletValveGeometry(pipenames(PIPENUM),4);

