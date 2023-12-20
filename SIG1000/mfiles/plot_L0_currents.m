clear all
close all
% stages of processing
% 1) define deployment number:
deploy  = 1
% 2) raw data input directory & filename convention:
rootDIR = ['/Users/derekgrimes/OneDriveUNCW/DATA/BOEM/FPSD1_Sig1000/mat_data/'];
fRoot   = ['S103071A010_FPS1_'];
% 3) output directory:
outRoot = ['/Users/derekgrimes/OneDriveUNCW/SIG1000/'];
% 4) output data file prefix:
filePrefix= sprintf('SIG_00103071_DEP%d_FPSC0_',deploy);
L0FRoot  = sprintf('%sL0_',filePrefix);
%
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
L1FRoot   = sprintf('%sL1',filePrefix);

load(



xm = 2;
ym = 2;
pw = 6;
ph = 3;
ag = 0.5;
cw = 0.25;
%
ppos1 = [xm ym pw ph];
cbpos1= [ppos1(1:2)+[ppos1(3) 0] cw 0.7*ph];
ppos2 = [xm ym+ph+ag pw ph];
cbpos2= [ppos2(1:2)+[ppos2(3) 0] cw 0.7*ph];
ppos3 = [xm ym+2*ph+2*ag pw ph];
cbpos3= [ppos3(1:2)+[ppos3(3) 0] cw 0.7*ph];
ppos4 = [xm ym+3*ph+3*ag pw ph];
cbpos4= [ppos4(1:2)+[ppos4(3) 0] cw 0.7*ph];
