function [KeyDownIdx,KeyUpIdx, DurAttacks, DurNotesInS, VelPeakIdxPos,VelPeakIdxNeg] = DetectVelocityPeaks_func(velo,tvec, FNAME)
pref = FNAME(1:3);
SR = 51.2e3;

Min = min(velo);
Max = max(velo);
if strcmp(pref,'A37')
    [VelPeakValPos, VelPeakIdxPos] = findpeaks(max(velo,0),'MinPeakProminence',0.087616 , 'MinPeakDistance', fix(51.2e3*1.000));    
else
    [VelPeakValPos, VelPeakIdxPos] = findpeaks(max(velo,0),'MinPeakProminence',0.5*Max , 'MinPeakDistance', fix(51.2e3*1.000));
end

[VelPeakValNeg, VelPeakIdxNeg] = findpeaks(max(-velo,0),'MinPeakProminence',-0.5*Min, 'MinPeakDistance', fix(51.2e3*1.000) );
VelPeakValNeg = -VelPeakValNeg;

% PARSE Upwards peaks in case of key bounce velocity
THR = fix(51.2e3*0.030);
for idx = 1:length(VelPeakIdxNeg)
    
    if ( VelPeakIdxPos(idx)- VelPeakIdxNeg(idx))<THR
        VelPeakIdxPos(idx) = [];
        VelPeakValPos(idx) = [];
    end 
end



if 0
figure(5); clf;
plot( (velo)  );
hold on;
plot(VelPeakIdxPos, VelPeakValPos,'ro');
plot(VelPeakIdxNeg, VelPeakValNeg,'r*');
end
% =========================================
if ne(length(VelPeakIdxNeg),length(VelPeakIdxPos))
    fprintf("NegPeaks number is not equal not PosPeak number\n")
    VelPeakValNeg = VelPeakValNeg(1:end-1);
    VelPeakIdxNeg = VelPeakIdxNeg(1:end-1);
end 
% mav = dsp.MovingAverage(fix(51.2e3*0.1e-3));

try
    DurNotesInS = tvec(VelPeakIdxPos-VelPeakIdxNeg);
catch
    DurNotesInS = NaN;
end

% DURATION OF KEY MOVING DOWNWARDS 

NP = length(VelPeakIdxNeg); % number of peaks

Tbefore = 0.200; % Assess noise starting at (s) before the lower peak
Tscan   = 0.100; % During this much time

TdownInit = zeros(length(VelPeakIdxNeg),1);
TdownEnd  = zeros(length(VelPeakIdxNeg),1);

% ASSUME the key will be moving downards less than this time:
Tmove = 0.100;
THRDWN = 0.05;
% % % figure(6); clf; %hold on;    
for idx = 1 : NP
    
% ASSES PRECEEDING NOISE 
seg = min(0, velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore)  :VelPeakIdxNeg(idx)-fix(SR*Tbefore)+fix(SR*Tscan) ) );
% AMPLITUDE:
PP = peak2peak( seg )*1;
STD = std(seg);

clf;
% % % % % % % plot(           velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore)  :VelPeakIdxNeg(idx)-fix(SR*Tbefore)+fix(SR*Tscan)   ));
% % plot(           velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore)  :VelPeakIdxNeg(idx)+fix(SR*0.060)   ));

% START NEGATIVE VELOCITY =============================

% InitDownIdx =  VelPeakIdxNeg(idx)-fix(SR*Tbefore) + find( (velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore): VelPeakIdxNeg(idx) )+PP)<0, 1, 'first');
% InitDownIdx =  VelPeakIdxNeg(idx)-fix(SR*Tbefore) + find( ( velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore): VelPeakIdxNeg(idx) ) + 20*PP )<0, 1, 'first');
InitDownIdx =  VelPeakIdxNeg(idx)-fix(SR*Tbefore) + find( ( min(0,velo( VelPeakIdxNeg(idx)-fix(SR*Tbefore): VelPeakIdxNeg(idx) )  ) )<THRDWN*(VelPeakValNeg(idx)), 1, 'first');
KeyDownIdx(idx) = InitDownIdx;
TdownInit(idx)  = tvec(InitDownIdx);


% END OF NEGATIVE VELOCITY  =============================

% segment = velo(VelPeakIdxNeg(idx) :VelPeakIdxNeg(idx)+fix(SR*Tmove) ); 
% EndDownInd    = VelPeakIdxNeg(idx) + find( velo(VelPeakIdxNeg(idx) :VelPeakIdxNeg(idx)+fix(SR*Tmove) ) <0, 1, 'last'  );

segment = velo( KeyDownIdx(idx) : KeyDownIdx(idx) + fix(SR*0.060) );
segment(1:fix(SR*0.013)) = [-1];
% FirstPosCrossingIdx = find( segment>0, 1, 'first');
FirstPosCrossingIdx = find( segment>THRDWN*(VelPeakValNeg(idx)), 1, 'first');
EndDownInd = KeyDownIdx(idx) + FirstPosCrossingIdx;
KeyUpIdx(idx) = EndDownInd;
TdownEnd(idx) = tvec(EndDownInd);

% % % hold on;
% % % plot([1,1]*(KeyDownIdx(idx)-VelPeakIdxNeg(idx)+ fix(SR*Tbefore) ),0.6*[-1,1],'--r');
% % % plot([1,1]*(KeyUpIdx(idx)-VelPeakIdxNeg(idx)+ fix(SR*Tbefore) ),0.6*[-1,1],'--k');
% % % hold off;
% % % pause();

end

% SO the key attack times are the differences of those two
DurAttacks = TdownEnd - TdownInit;



