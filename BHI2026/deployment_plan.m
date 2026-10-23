%
%
clear all
close all


%% Pre-deployment plan
%% ECO ADCP's only...
UTM_E_ADCP = [226180  225212];
UTM_N_ADCP = [3749481 3748954];

%% Pressure sensors only...
UTM_E = [226180  226180 225712 224712];
UTM_N = [3748981 3749981 3748754 3749154];
% $$$ 
proj = projcrs(32618);
[Lat,Lon]           = projinv(proj,UTM_E,UTM_N);
[Lat_ADCP,Lon_ADCP] = projinv(proj,UTM_E_ADCP,UTM_N_ADCP);

figure,
geoplot(Lat,Lon,'*'), hold on,
geoplot(Lat_ADCP,Lon_ADCP,'*'),
geobasemap('satellite')

return

fin  = '/Users/derekgrimes/git/FPS_PhysOcn/BHI2026/deploy_locations.csv';
fid  = fopen(fin);
format = ['%d%f%f%f%s%*s%*s%s', repmat('%*s',1,20)]
data = textscan(fid,format,'HeaderLines',1,'delimiter',',')
Name= data{6};
Lat = data{2};
Lon = data{3};
% $$$ Time= datetime(char(data{5}))
figure,
geoplot(Lat,Lon,'*'),
geobasemap('satellite')