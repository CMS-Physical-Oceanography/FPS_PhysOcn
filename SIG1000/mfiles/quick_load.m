% get a quick look at deployment data
%
rootDIR = '/Users/derekgrimes/OneDriveUNCW/DATA/BOEM/FPSD1_Sig1000/mat_data/';
files   = dir([rootDIR,'*.mat']);
%
load([rootDIR,filesep,files(1).name],'Data','Config')
%
% extract all variables starting with "burst"
vars  = fields(Data);
burst = regexp(vars,'(?<=Burst_).*','match');
echo1 = regexp(vars,'(?<=Echo1Bin1_1000kHz_).*','match');
echo2 = regexp(vars,'(?<=Echo2Bin1_1000kHz_).*','match');
for ii=1:length(burst);
    varBurst = char(burst{ii});
    varEcho1 = ['Echo1_',char(echo1{ii})];
    varEcho2 = ['Echo2_',char(echo2{ii})];
    if ~isempty(varBurst)
        eval([varBurst,' = Data.',char(vars{ii}),';'])
    end

    if ~isempty(varEcho1)
        eval([varEcho1,' = Data.',char(vars{ii}),';'])
    end

    if ~isempty(varEcho2)
        eval([varEcho2,' = Data.',char(vars{ii}),';'])
    end
end
%
% get bin depths
z0  = Config.Burst_BlankingDistance;
dz  = Config.Burst_CellSize;
nz  = double(Config.Burst_NCells);
z   = z0+[0:(nz-1)]*dz+0.5*dz;
%
% find bins that are in the water
inW = (z<(Pressure-dz));
%
dt = (Time(2)-Time(1))*86400;
%
% 5 minute avg;
N  = round(300/dt);
w  = hanning(N); w = w./sum(w);
%
inW = conv2(w,1,double(inW),'same');
inW = inW>0.5;
%
vE  = squeeze(Velocity_ENU(:,1,:));
vN  = squeeze(Velocity_ENU(:,2,:));
vU  = squeeze(Velocity_ENU(:,3,:));
err = squeeze(Velocity_ENU(:,4,:));
%
%
vEavg = conv2(w,1,vE,'same'); vEavg(~inW)=nan;
vNavg = conv2(w,1,vN,'same'); vNavg(~inW)=nan;
vUavg = conv2(w,1,vU,'same'); vUavg(~inW)=nan;
erravg = conv2(w,1,err,'same'); erravg(~inW)=nan;
figure, imagesc(Time,z,vEavg'),title('East'),colormap(cmocean('balance')),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,caxis([-0.5 0.5])
figure, imagesc(Time,z,vNavg'),title('North'),colormap(cmocean('balance')),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,caxis([-0.25 0.25])
figure, imagesc(Time,z,vUavg'),title('Up'),colormap(cmocean('balance')),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,caxis([-0.05 0.05])
figure, imagesc(Time,z,erravg'),title('Error'),colormap(cmocean('balance')),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,caxis([-0.01 0.01])
%
%
% get bin depths
z0E  = Config.EchoSounder_BlankingDistance;
dzE  = Config.EchoSounder_CellSize;
nzE  = double(Config.EchoSounder_NCells);
zE   = z0E+[0:(nzE-1)]*dzE+0.5*dzE;
%
E1  = conv2(w,1,Echo1_Echo);
E2  = conv2(w,1,Echo2_Echo);
figure, imagesc(Echo1_Time,zE,E1'),title('Echo1'),colormap(bone),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,

figure, imagesc(Echo2_Time,zE,E2'),title('Echo2'),colormap(bone),set(gca,'ydir','normal'),datetick('x','keepticks','keeplimits'),colorbar,
