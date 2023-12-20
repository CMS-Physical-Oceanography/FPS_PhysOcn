clear all
close all
% load input averaged file
fin = '/Users/derekgrimes/OneDriveUNCW/SIG1000/SIG_00103071_DEP1_FPSC0_L1.mat';
load(fin)
fout= '/Users/derekgrimes/PHY477/HW/HW7/FPS_UVP.csv';
%
% mask out flagged data
u = VelEast; u(~qcFlag)=nan;
v = VelNorth;v(~qcFlag)=nan;
Time = cellstr(datestr(Time));
EastVelocity = nanmean(u,2);
NorthVelocity= nanmean(v,2);
T   = table(Time,EastVelocity,NorthVelocity,Pressure);
writetable(T,fout)
% $$$ fid = fopen(fout,'w');
% $$$ fprintf(fid,'%s %f %f %f',M);
% $$$ 
% $$$ for j = 1:length(Time)
% $$$     fprintf(fid,'%s %f %f %f',datestr(time(j)),nanmean(u(j,:)),nanmean(v(j,:)),Pressure(j));
% $$$ end