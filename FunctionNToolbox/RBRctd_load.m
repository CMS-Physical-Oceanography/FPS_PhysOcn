function [RBRctd, filename] = RBRctd_load(file_path,save_to_path,depnum,stationID,stime,etime);
% [RBRctd,filename] = RBRctd_load(file_path,save_to_path,depnum,stationID,stime,etime);
%
%=========================================================================
% RBRctd_load reads the .rsk file output from the Ruskin.exe
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
% conductivity, salinity, temperature and pressure, and their units. 
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
RBRctd.conductivity=RSK.data.values(:,1);
RBRctd.temperature=RSK.data.values(:,2);
RBRctd.pressure=RSK.data.values(:,3);
RBRctd.time=RSK.data.tstamp;
RBRctd.SN=string(rsk.instruments.serialID);
RBRctd.units= {'temp = deg C'; 'pres = dbar'};
RBRctd.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};

sal = real(sw_salt(RBRctd.conductivity./sw_c3515,RBRctd.temperature,RBRctd.pressure));
%%
% Trim the beginning and ends of the files:
aa = find(RBRctd.time>=stime & RBRctd.time <=etime);
time = RBRctd.time(aa);
temperature = RBRctd.temperature(aa);
pressure = RBRctd.pressure(aa);
conductivity = RBRctd.conductivity(aa);
salinity = sal(aa);
units = RBRctd.units;
SN = RBRctd.SN;
notes = RBRctd.notes;

%% Raw figures
  figure(); clf;
     subplot(311);
        plot(time,pressure, 'LineWidth',1.5,'Color','blue');
        title(sprintf('RBR-CTD data %s',RBRctd.SN)); ylabel('dbar'); grid on;
        set(gca,'XTickLabel',[]); axis tight; 
    subplot(312);
        plot(time,temperature, 'LineWidth',1.5,'Color','red');
        ylabel('^oC'); grid on; axis tight; datetick('x','mm/dd','keepticks', 'keeplimits');
    subplot(313);
        plot(time,salinity, 'LineWidth',1.5,'Color','black');
        ylabel('^oC'); grid on; axis tight; datetick('x','mm/dd','keepticks', 'keeplimits');  
     
%% Saving the file name
filename=[save_to_path,'RBR_',sprintf('%08d',str2num(RBRctd.SN)),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
save(filename,'SN','time','temperature','pressure','salinity','conductivity','units','notes'); 
end
end









