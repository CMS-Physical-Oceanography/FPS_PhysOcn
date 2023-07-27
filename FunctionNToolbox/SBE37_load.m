function [SBE37] = SBE37_load(deploypath)

%=========================================================================
% This functions inputs the file path of the SeaTerm .asv output file and 
% outputs a Data structure with temperature, conductivity, pressure, &
% time into L0processing folder. 
%=========================================================================
sbe37Folder=dir([deploypath, filesep, 'SBE*']);
inpath=[deploypath, filesep, sbe37Folder(1).name, filesep];

source=fullfile(inpath,'*.asc');
filedir=dir(source);

file_path=fullfile(inpath,filedir(1).name);
fileID=fopen(file_path,'r');
formatspec='%f %f %f %f %s %s';

numLinesToSkip=8;

for i=1:numLinesToSkip
    dum=fgetl(fileID);
end
%Pulling the SN from the .asc file
SNline=fgetl(fileID);
SN=regexp(SNline,'(?<=NO.\s).*(?=\s30)','match');
SNstr = strtrim(SN);

%Skip lines down to data listed for the deployment
numLine= 49;
for i=1:numLine
    dum2=fgetl(fileID);
end

tline=fgetl(fileID);
sample_int=str2double(regexp(tline,'(?<==\s).*(?=\sseconds)','match'));
tline2=fgetl(fileID);

dataArray = textscan(fileID, formatspec , 'Delimiter', ',');
date_time=[char(dataArray{5}), repmat(' ',size(dataArray{1},1)), char(dataArray{6})];
time=datenum(date_time);
temp=dataArray{1,1}; %deg C
conductivity=dataArray{1,2}; %Siemens/m (S/m)
pres=dataArray{1,3};
sal=dataArray{1,4};

SBE37.temperature=temp;
SBE37.conductivity=conductivity;
SBE37.pressure=pres;
SBE37.salinity=sal;
SBE37.time=time;
SBE37.SN=SNstr;
SBE37.units={'temp = deg C'; 'cond = S/m'; 'pres = dB'; 'sal = PSU'};
%% Raw Data Figure
figure(); clf;
subplot(411)
    plot(time, temp,'LineWidth',1.5,'Color','blue'); ylabel('Temp (Deg C)');
    title('Raw SBE37 Data'); grid on; set(gca,'XTickLabel',[]); axis tight;
subplot(412)
    plot(time, conductivity, 'LineWidth', 1.5, 'Color', 'red'); ylabel('Cond (S/m)');
    grid on; axis tight; set(gca,'XTickLabel',[]);
subplot(413)
    plot(time, pres, 'LineWidth',1.5, 'Color', 'green'); ylabel('Pres (dB)'); grid on;
    axis tight; set(gca,'XTickLabel',[]);
subplot(414)
    plot(time, sal, 'LineWidth',1.5, 'Color', 'black'); ylabel('Sal (ppt)'); grid on;
    axis tight; datetick('x','dd/mm', 'keepticks', 'keeplimits');
    
%% Saving the file name
% Creating L0 Folder within the BOEMTest and Crate Folder
L0Folder=mkdir(deploypath, 'L0processing')

outpath=[deploypath, '..']; cd(outpath); outpath=pwd;
date=datestr(datetime(time(1), 'Convertfrom','datenum'),'yyyy_mm_dd');
%Linking the file to be saved in L0proc under BOEMtest folder
L0Folder=dir([outpath, filesep, 'L0*']);
filename=[outpath, filesep, L0Folder(1).name, filesep,'SBE37_',char(SBE37.SN),'_', date,'_L0.mat'];
save(filename, 'SBE37'); 



