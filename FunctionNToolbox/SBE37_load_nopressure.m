function [SBE37,filename] = SBE37_load_nopressure(file_path,save_to_path,depnum,stationID,stime,etime)
%[SBE37,filename] = SBE37_load_nopressure(file_path,save_to_path,depnum,stationID,stime,etime)
%
%=========================================================================
% This functions inputs the file path of the SeaTerm .asc output file and 
% outputs a Data structure with temperature, conductivity, pressure, &
% time into L0processing folder. 
%
% NOTE: This one does not include the pressure measurements
%=========================================================================
isthere=dir(file_path);
if isempty(isthere) == 1
disp(strcat('Can''t find  ',file_path))
else
    
fileID=fopen(file_path,'r');
% Really annoying - changes whether or not sound speed is calculated as
% well as the indexing below
formatspec='%f %f %f %f %s %s';
% 
% Read the header and stop at line that starts with '*END*'
str='*';
SN_pattern = 'SERIAL NO\. (\d+)';
while (~strncmp(str,'*END*',5));
   str=fgetl(fileID);
   if isempty(regexp(str, SN_pattern, 'tokens', 'once')) == 1
   else
       SN = cell2mat(regexp(str, SN_pattern, 'tokens', 'once'));
   end
end

%Then there are 3 lines to skip:
numLinesToSkip=3;
for i=1:numLinesToSkip
    dum=fgetl(fileID);
end

dataArray = textscan(fileID, formatspec , 'Delimiter', ',');
fclose(fileID)

date_time = char(strcat(dataArray{5},{' '},dataArray{6}));

time=datenum(date_time);
temp=dataArray{1,1}; %deg C
conductivity=dataArray{1,2}; %Siemens/m (S/m)
pres=NaN(size(conductivity));
sal=dataArray{1,3};

SBE37.temperature=temp;
SBE37.conductivity=conductivity;
SBE37.pressure=pres;
SBE37.salinity=sal;
SBE37.time=time;
SBE37.SN=SN;
SBE37.units={'temp = deg C'; 'cond = S/m'; 'pres = dB'; 'sal = PSU'};
SBE37.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
%%
% Trim the beginning and ends of the files:
aa = find(SBE37.time>=stime & SBE37.time <=etime);
time = SBE37.time(aa);
temperature = SBE37.temperature(aa);
conductivity = SBE37.conductivity(aa);
pressure = SBE37.pressure(aa);
salinity = SBE37.salinity(aa);
units = SBE37.units;
SN = SBE37.SN;
notes = SBE37.notes;

%% Raw Data Figure
figure(); clf;
subplot(411)
    plot(time, temperature,'LineWidth',1.5,'Color','blue'); ylabel('Temp (Deg C)');
    title(sprintf('Raw SBE37 Data %s', SBE37.SN)); grid on; set(gca,'XTickLabel',[]); axis tight; 
subplot(412)
    plot(time, conductivity, 'LineWidth', 1.5, 'Color', 'red'); ylabel('Cond (S/m)');
    grid on; axis tight; set(gca,'XTickLabel',[]); ylim([0 7.5])
subplot(413)
    plot(time, pressure, 'LineWidth',1.5, 'Color', 'green'); ylabel('Pres (dB)'); grid on;
    axis tight; set(gca,'XTickLabel',[]);
subplot(414)
    plot(time, salinity, 'LineWidth',1.5, 'Color', 'black'); ylabel('Sal (PSU)'); grid on;
    axis tight; datetick('x','mm/dd', 'keepticks', 'keeplimits');
    
%% Saving the file name
filename=[save_to_path,'SBE_',sprintf('%08d',str2num(SN)),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
save(filename,'SN','time','temperature','conductivity','salinity','pressure','units','notes'); 
end



