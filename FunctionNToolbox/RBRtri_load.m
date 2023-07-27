function [RBRtri] = RBRtri_load(deploypath)

%=========================================================================
% RBRtri_load reads the .rsk file output from the Ruskin.exe
% software. The inputs to the function are the full filename ".rsk" and the 
% file directory to the RSKtools. 

%|RSKtools| is RBR's open source Matlab toolbox for reading,
% visualizing, and post-processing RBR logger data. It provides
% high-speed access to large RSK data files. Users may plot data as a
% time series or as depth profiles using tailored plotting
% utilities. If the tool box is not already downloaded, it
%can be found at this link: https://rbr-global.com/support/matlab-tools/
%
% The output for this file is a data structure that contains time,
% chlorophyll-a, FDOM, and turbidity data and units. 
%
%If a user needs more clarification on certain parts of the tool box, type
%"help RSKtools" in the command window. 
%=========================================================================

%Creating path to RBRtri within BOEMtest folder
rbrtriFolder=dir([deploypath, filesep, 'RBRtri*']);
inpath=[deploypath, filesep, rbrtriFolder(1).name];
source=fullfile(inpath,'*.rsk');
filedir=dir(source);
file_path=fullfile(inpath,filedir(1).name);

%Reading raw .rsk file into matlab structure rsk
rsk=RSKopen(file_path);

% Specify the time frame to read data from the RBR
RSK=RSKreaddata(rsk);
%rsk = RSKreaddata(rsk, 't1', t1, 't2', t2);

%Data Output Structure
RBRtri.chlorophyll_a=RSK.data.values(:,1);
RBRtri.FDOM=RSK.data.values(:,2);
RBRtri.turbidity=RSK.data.values(:,3);
RBRtri.time=RSK.data.tstamp;
RBRtri.SN=num2str(rsk.instruments.serialID);
RBRtri.units= {'chlorophyll-a = ug/l'; 'FDOM = ppb'; 'turb = FTU'};
%% Raw figures
  figure(); clf;
     subplot(311);
        plot(RBRtri.time,RBRtri.chlorophyll_a, 'LineWidth',1.5,'Color','blue');
        title('RBR Tridente Raw Data'); ylabel('phyll-a (ug/l)'); grid on;
        set(gca,'XTickLabel',[]); axis tight
    subplot(312);
        plot(RBRtri.time, RBRtri.FDOM, 'LineWidth',1.5,'Color','green');
        ylabel('FDOM (ppb)'); grid on; set(gca,'XTickLabel',[]); axis tight
    subplot(313);
        plot(RBRtri.time, RBRtri.turbidity, 'LineWidth',1.5,'Color','red');
        ylabel('turb (FTU)'); grid on; axis tight; datetick('x','dd/mm','keepticks', 'keeplimits')
     
%% %Saving the file to outpath directory with file name and start date of
%deployment

% Creating L0 Folder within the BOEMTest and Crate Folder
L0Folder=mkdir(deploypath, 'L0processing')

date=datestr(RBRtri.time(1),'yyyy_mm_dd');
%Linking the file to be saved in L0proc under BOEMtest
outpath=[inpath, filesep, '..'];
cd(outpath); outpath=pwd;
L0Folder=dir([outpath, filesep, 'L0*']);
filename=[outpath, filesep, L0Folder(1).name, filesep, 'RBRtri_', RBRtri.SN,'_', date,'_L0.mat'];
save(filename, 'RBRtri'); 
end









