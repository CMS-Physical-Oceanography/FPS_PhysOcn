function [RBRduet, filename] = RBRduet_load(file_path,save_to_path,depnum,stationID,stime,etime);
% [RBRduet,filename] = RBRduet_load(file_path,save_to_path,depnum,stationID,stime,etime);
%
%=========================================================================
% RBRduet_load reads the .rsk file output from the Ruskin.exe
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
% temperature and pressure, and their units. 
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
RBRduet.temperature=RSK.data.values(:,1);
RBRduet.pressure=RSK.data.values(:,2);
RBRduet.time=RSK.data.tstamp;
RBRduet.SN=string(rsk.instruments.serialID);
RBRduet.units= {'temp = deg C'; 'pres = dbar'};
RBRduet.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
%%
% Trim the beginning and ends of the files:
aa = find(RBRduet.time>=stime & RBRduet.time <=etime);
time = RBRduet.time(aa);
temperature = RBRduet.temperature(aa);
pressure = RBRduet.pressure(aa);
units = RBRduet.units;
SN = RBRduet.SN;
notes = RBRduet.notes;

%% Raw figures
  figure(); clf;
     subplot(211);
        plot(time,pressure, 'LineWidth',1.5,'Color','blue');
        title(sprintf('RBR-duet data %s',RBRduet.SN)); ylabel('dbar'); grid on;
        set(gca,'XTickLabel',[]); axis tight; 
    subplot(212);
        plot(time,temperature, 'LineWidth',1.5,'Color','red');
        ylabel('^oC'); grid on; axis tight; datetick('x','mm/dd','keepticks', 'keeplimits');
        
     
%% Saving the file name
filename=[save_to_path,'RBR_',sprintf('%08d',str2num(RBRduet.SN)),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
save(filename,'SN','time','temperature','pressure','units','notes'); 
end
end









