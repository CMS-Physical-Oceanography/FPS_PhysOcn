function [RBRsoloT] = RBRsoloT_load(deploypath)
% RBRsoloT_load loads in the path (deploypath) to the .rsk file for the
% RBRsoloT and makes a raw plot of temperature data. 

%Loading in the .rsk file from RBRsoloT folder
RBRtempFolder=dir([deploypath, filesep, 'RBRs*']);
inpath=[deploypath, filesep, RBRtempFolder.name];
soloTfile=dir(fullfile(inpath, '*.rsk'));
file_path=fullfile(inpath, soloTfile.name);

%RSKtools toolbox used to process data into matlab from the .rsk file
rsk=RSKopen(file_path);

%Finding the start and end time of the data within the file
startTime=rsk.epochs.startTime;
endTime=rsk.epochs.endTime;

% Specify the time frame to read data from the RBR
%the # after startTime regards how many days the function RSKread
RSK=RSKreaddata(rsk);

%Data Output Structure
RBRsoloT.SN=num2str(RSK.instruments.serialID);
RBRsoloT.temp=RSK.data.values;
RBRsoloT.time=RSK.data.tstamp;
RBRsoloT.units='deg C';

%% Raw figures
figure;
plot(RBRsoloT.time, RBRsoloT.temp, 'LineWidth',1.5,'Color', 'red');
title('Raw RBRsoloT Temp'); ylabel('Deg C'); grid on; axis tight;
datetick('x','dd/mm','keepticks', 'keeplimits')
%% Saving file to L0 proc folder within the Correct BOEMTest and Crate

% Creating L0 Folder within the BOEMTest and Crate Folder
L0Folder=mkdir(deploypath, 'L0processing')

%create a date time of the first time stamp for logging reference
date=datestr(RBRsoloT.time(1),'yyyy_mm_dd');
%Linking the file to be saved in L0proc under BOEMtest
outpath=[inpath, filesep, '..'];
cd(outpath); outpath=pwd;
L0Folder=dir([outpath, filesep, 'L0*']);
filename=[outpath, filesep, L0Folder(1).name, filesep, 'RBRsoloT_', RBRsoloT.SN ,'_',date ,'_L0.mat'];
save(filename, 'RBRsoloT');









end