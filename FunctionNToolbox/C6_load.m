function [C6,filename] = C6_load(file_path,save_to_path,depnum,stationID,stime,etime)
% [C6,filename] = C6_load(file_path,save_to_path,depnum,stationID,stime,etime)
%
% ========================================================================%
%C6_load takes the raw data .csv file from the C6 fluoremeter and inputs it
%into a structure called C6. 
% Function Input:
%   file_path- full path to the deployment folder with .csv filename. 
%   save_to_path- full path to the folder to save the L0 file (NONAME). 
%   depnum - 1 through 8
%   stationID - 'FPSC', 'FPSE1', 'FPSS1', etc.
%
% Function Output:
%   C6 - data struture containing the instruments SN, CDOM, Chlorophyll a,
%        Turbidity, Depth, and Temp.
% ========================================================================%

isthere=dir(file_path);
if isempty(isthere) == 1
disp(strcat('Can''t find  ',file_path))
else
    
%Using fopen command to read line by line of the .csv file in notepad++
fileID=fopen(file_path,'r');
%Pulling the Serial Number from the .csv txt file
SN_line= fgetl(fileID);
SN = regexp(SN_line, '(\d+)', 'match');

 numLinestoSkip = 8;
 for i= 1:numLinestoSkip
     dum=fgetl(fileID);
 end
%Specify text labels to implement into the textscan function to read the data lines
% of the C6 .csv txt file
formatspec = '%s %f %f %f %f %f';
dataArray = textscan(fileID, formatspec , 'Delimiter', ',');
time=datenum(dataArray{1,1});
CDOM = dataArray{1,2};
chlorophyll_a = dataArray{1,3};
turbidity = dataArray{1,4};
depth = dataArray{1,5};
temperature = dataArray{1,6};
units = {'CDOM = ug/L' ; 'Chlorophyll-a = ug/L' ; 'Turbidity = NTU'; 'Depth = meters'; 'Temperature = Deg C'};
fclose(fileID)

%Create the Data structure
C6.time =time;
C6.CDOM = CDOM;
C6.chlorophyll_a = chlorophyll_a;
C6.turbidity = turbidity;
C6.depth = depth;
C6.temperature = temperature;
C6.units = units; 
C6.SN = cell2mat(SN);
C6.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
%%
% Trim the beginning and ends of the files:
aa = find(C6.time>=stime & C6.time <=etime);
time = C6.time(aa);
CDOM = C6.CDOM(aa);
chlorophyll_a = C6.chlorophyll_a(aa);
turbidity = C6.turbidity(aa);
depth = C6.depth(aa);
temperature = C6.temperature(aa);
units = C6.units;
SN = C6.SN;
notes = C6.notes;

% Raw data plots

figure
subplot(511)
    plot(time, CDOM, 'b','LineWidth',1.5); grid on; axis tight;
    ylabel('CDOM (ug/L)'); set(gca,'XtickLabels', []); title(sprintf('C6 Raw Data: %s', C6.SN));
    ylim([0 25]);
subplot(512)
    plot(time, chlorophyll_a, 'g','LineWidth',1.5); grid on; axis tight;
    ylabel('Chlorophyll-a (ug/L)'); ylim([0 25]); set(gca,'XTickLabel',[]);
subplot(513)
    plot(time, turbidity, 'k', 'LineWidth',1.5); axis tight; grid on;
    ylabel('Turbidity (NTU)'); set(gca, 'XTickLabel', []);
subplot(514)
    plot(time, depth, 'm','LineWidth',1.5); grid on; axis tight;
    ylabel('Depth (m)'); set(gca, 'XTickLabel', []);
subplot(515)
    plot(time, temperature, 'r','LineWidth',1.5); axis tight; grid on
    datetick('x','mm/dd','keeplimits','keepticks'); ylabel('Temperature (^oC)');

%% Saving the file name
filename=[save_to_path,'0C6_',char(C6.SN),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
save(filename,'SN','time','CDOM','chlorophyll_a','turbidity','depth','temperature','units','notes'); 
end


