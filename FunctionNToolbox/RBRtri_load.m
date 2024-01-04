function [RBRtri, filename] = RBRtri_load(file_path,save_to_path,depnum,stationID,stime,etime);
% [RBRtri, filename] = RBRtri_load(file_path,save_to_path,depnum,stationID,stime,etime);
%
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

isthere=dir(file_path);
if isempty(isthere) == 1
disp(strcat('Can''t find  ',file_path))
else
    
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
RBRtri.SN=string(rsk.instruments.serialID);
RBRtri.units= {'chlorophyll-a = ug/l'; 'FDOM = ppb'; 'turb = FTU'};
RBRtri.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
%%
% Trim the beginning and ends of the files:
aa = find(RBRtri.time>=stime & RBRtri.time <=etime);
time = RBRtri.time(aa);
chlorophyll_a = RBRtri.chlorophyll_a(aa);
FDOM = RBRtri.FDOM(aa);
turbidity = RBRtri.turbidity(aa);
units = RBRtri.units;
SN = RBRtri.SN;
notes = RBRtri.notes;

%% Raw figures
  figure(); clf;
     subplot(311);
        plot(time,chlorophyll_a, 'LineWidth',1.5,'Color','blue');
        title(sprintf('RBR Tridente Raw Data %s',RBRtri.SN)); ylabel('phyll-a (ug/l)'); grid on;
        set(gca,'XTickLabel',[]); axis tight; 
    subplot(312);
        plot(time, FDOM, 'LineWidth',1.5,'Color','green');
        ylabel('FDOM (ppb)'); grid on; set(gca,'XTickLabel',[]); axis tight; 
    subplot(313);
        plot(time, turbidity, 'LineWidth',1.5,'Color','red');
        ylabel('turb (FTU)'); grid on; axis tight; datetick('x','dd/mm','keepticks', 'keeplimits');
        
     
%% %Saving the file to outpath directory with file name and start date of
filename=[save_to_path,'RBR_',sprintf('%08d',str2num(RBRtri.SN)),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
save(filename,'SN','time','chlorophyll_a','FDOM','turbidity','units','notes'); 
end
end









