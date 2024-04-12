Config.ATM_Time = ATM_Time;
Config.ATM_Pressure=ATM_Pressure;
%
% average the start/end pressures
ATM_Pressure = mean(ATM_Pressure);
%
% get bin sizing etc for currents
blnk = Config.Burst_BlankingDistance;
binw = Config.Burst_CellSize;
Nc   = double(Config.Burst_NCells);
mab  = hab+blnk+binw.*[1:Nc];
%
%  get bin sizing etc for echo
blnkE= Config.EchoSounder_BlankingDistance;
binwE= Config.EchoSounder_CellSize;
NcE  = double(Config.EchoSounder_NCells);
mabE = hab+blnkE+binwE.*[1:NcE];
%
% save the config info
outFile = [outRoot, outFRoot, 'config.mat'];
save(outFile,'Config','Descriptions')
%
% limit archive file to 24hr length
maxDuration = 24*3600;
%
% get sample rate
fs = double(Config.Burst_SamplingRate);
%
% maximum number of samples per file
Nmax = fs*maxDuration;
% number of measurements in current file
N    = 0;
%
% flag to load next file
loadFlag = 1;
ii = 0;
outNf = 0;
fprintf('\n \n')
while ii<=Nf
    %
    if loadFlag
        ii  = ii+1;
        loadFlag=0;
        fileName = [fRoot,num2str(ii),'.mat'];
        fin = [files(ii).folder,'/',fileName];
        fprintf(['pre-processing:       %s \n'], fileName)
        load(fin)
        % use 5th beam time; the other beams are offset by dt = 1/(2*fs)
        t    = Data.IBurst_Time+time_shift;
        is   = find(t>=deployTime,1,'first');
        nt   = length(t);
        % check if instrument was pulled
        iStop= find(t>=recoverTime,1,'first');
        if ~isempty(iStop)
            disp(['trimming data at recovery time: nt = ',num2str(nt),'--> ',num2str(iStop-1)]);
            nt = iStop-1;
        end
    end
    %
    if N==0
        % initialize output structure
        out  = struct('Time',[],'Heading',[],'Pitch',[],'Roll',[],'Temperature',[],'Pressure',[],'Battery',[],...
                      'VelBeam1',[],'VelBeam2',[],'VelBeam3',[],'VelBeam4',[],'VelBeam5',[],...
                      'AmpBeam1',[],'AmpBeam2',[],'AmpBeam3',[],'AmpBeam4',[],'AmpBeam5',[],...
                      'CorBeam1',[],'CorBeam2',[],'CorBeam3',[],'CorBeam4',[],'CorBeam5',[],...,
                      'VelEast',[],'VelNorth',[],'VelUp',[],'VelError',[],...
                      'mab',mab, 'mab_Echo',mabE,'Echo1',[],'Echo2',[]);
    end
    %
    if N+nt-(is-1)<=Nmax
        ie   = nt;
        loadFlag = 1;
    else
        ie   = Nmax-N-(is-1);
        loadFlag = 0;
    end
    %
    out.Time       (N+1:N+(ie-(is-1)),1)     = Data.IBurst_Time(is:ie,1)+time_shift;
    out.Heading    (N+1:N+(ie-(is-1)),1)     = Data.Burst_Heading(is:ie,1);
    out.Pitch      (N+1:N+(ie-(is-1)),1)     = Data.Burst_Pitch(is:ie,1);
    out.Roll       (N+1:N+(ie-(is-1)),1)     = Data.Burst_Roll(is:ie,1);
    out.Temperature(N+1:N+(ie-(is-1)),1:2)  = [Data.IBurst_Temperature(is:ie,1), Data.Burst_Temperature(is:ie,1)];
    out.Pressure   (N+1:N+(ie-(is-1)),1:2)  = [Data.IBurst_Pressure(is:ie,1),    Data.Burst_Pressure(is:ie,1)] - ATM_Pressure;
    out.Battery    (N+1:N+(ie-(is-1)),1)     = Data.IBurst_Battery(is:ie,1);
    out.VelBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_VelBeam1(is:ie,1:Nc);
    out.VelBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_VelBeam2(is:ie,1:Nc);
    out.VelBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_VelBeam3(is:ie,1:Nc);
    out.VelBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_VelBeam4(is:ie,1:Nc);    
    out.VelBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_VelBeam5(is:ie,1:Nc);
    out.AmpBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_AmpBeam1(is:ie,1:Nc);
    out.AmpBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_AmpBeam2(is:ie,1:Nc);
    out.AmpBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_AmpBeam3(is:ie,1:Nc);
    out.AmpBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_AmpBeam4(is:ie,1:Nc);    
    out.AmpBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_AmpBeam5(is:ie,1:Nc);
    out.CorBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_CorBeam1(is:ie,1:Nc);
    out.CorBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_CorBeam2(is:ie,1:Nc);
    out.CorBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_CorBeam3(is:ie,1:Nc);
    out.CorBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_CorBeam4(is:ie,1:Nc);    
    out.CorBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_CorBeam5(is:ie,1:Nc);
    out.VelEast    (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_ENU(is:ie,1,1:Nc);
    out.VelNorth   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_ENU(is:ie,2,1:Nc);
    out.VelUp      (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_ENU(is:ie,3,1:Nc);
    out.VelError   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_ENU(is:ie,4,1:Nc);    
    %
    N = N+ie-(is-1);
    is= ie+1;
    %
    if N==Nmax || (ii==Nf & loadFlag)
        % process data
        minC  = min( cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        minA  = min( cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        maxRNG= out.Pressure(:,2)*cosd(25)-binw;
        qcFlag= (minC>50 & minA>30 & maxRNG>mab);
        %
        out.minCor = minC;
        out.minAmp = minA;
        out.maxRNG = maxRNG;
        out.qcFlag = qcFlag;
        %
        disp('not applying heading offset correction, velocities are relative to magnetic north/south')
% $$$         % now rotate to xyz and ENU
% $$$         [out,Config2,beam2xyz] = beam2earth_sig1000_DG_FRFcoords(out,Config,'');        
% $$$         out.HeadingOffset = HeadingOffset;
% $$$         out.beam2xyz = beam2xyz;
        %
        % save output
        outNf = outNf+1;
        outFileName = sprintf([outFRoot,'%03d.mat'],outNf);
        fout  = [outRoot,outFileName];
        fprintf(['saving output file:   %s \n'],outFileName)
        save(fout,'-struct','out')
        %
        N = 0;
        %
        if (ii==Nf & loadFlag )
        disp('done!')
        break
        end
    end
    %
    %
end


