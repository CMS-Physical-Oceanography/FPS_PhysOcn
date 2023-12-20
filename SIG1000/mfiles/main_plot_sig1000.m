clear all
close all
% stages of processing
% 1) define deployment number:
deploy    = 1;
filePrefix= sprintf('SIG_00103071_DEP%d_FPSC0_',deploy);
% 3) output directory:
outRoot   = ['/Users/derekgrimes/OneDriveUNCW/SIG1000/'];
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
L1FRoot   = sprintf('%sL1',filePrefix);
% 5) input time-averaged file
avgFile   = [outRoot,filesep,L1FRoot,'.mat'];
% 6) output tidal constituents file
tideFile  = [outRoot,filesep,L1FRoot,'_tidal_analysis.out'];
%
load(avgFile)
%
% estimate tidal components, sample rate in decimal hours
[Nt,Nb]  = size(VelEast);       
dt       = dtAvg/3600;
b2u      = (sum(~qcFlag)==0);
TideEast = nan(Nt,Nb);
TideNorth= nan(Nt,Nb);
fid = fopen(tideFile,'w');
fseek(fid,0,1);
fprintf(fid,'Tidal Analysis for deployment #%d pressure record:\n',deploy);
fclose(fid)
[~,~,~,TidePres] = t_tide(Pressure-mean(Pressure),'interval',dt,'output',tideFile);
for bin  = find(b2u);
    uv = VelEast(:,bin) + sqrt(-1).*VelNorth(:,bin);
    fid = fopen(tideFile,'a');
    fprintf(fid,['\n\n--------------------------------------\n',...
                 'velocity record in bin number: %d\n',...
                 '--------------------------------------\n'],bin);
    fclose(fid)
    [~,~,~,uv_tide] = t_tide(uv,'interval',dt,'output',tideFile);
    TideEast (:,bin) = real(uv_tide);
    TideNorth(:,bin) = imag(uv_tide);
end
%
% interpolate outside of "good" bins
indI = logical(ones(Nt,1)*b2u);
indJ = qcFlag & ~indI;
%
tt   = Time*ones(1,Nb);
zz   = ones(Nt,1)*mab;
F = scatteredInterpolant(tt(indI),zz(indI),TideEast(indI),'linear','nearest');
TideEast(indJ) = F(tt(indJ),zz(indJ));
F.Values = TideNorth(indI);
TideNorth(indJ) = F(tt(indJ),zz(indJ));
% TideEast(indJ) = interp2(tt(indI),zz(indI),TideEast(indI),tt(indJ),zz(indJ));