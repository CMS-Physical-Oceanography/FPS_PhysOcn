function [SB37] = SB37_process_deployment(path,filename)

%=========================================================================
% This functions inputs the file path of the SeaTerm .asv output file and 
% outputs a Data structure with temperature, conductivity, pressure, &
% time. 
%=========================================================================

filedir=append(path,'\',filename);
fileID=fopen(filedir,'r')
formatspec='%f %f %f %f %s %s';

numLinesToSkip=57;

for i=1:numLinesToSkip
    dum=fgetl(fileID);
end

tline=fgetl(fileID);
tline2=fgetl(fileID);
sample_int=str2double(regexp(tline2,'(?<==\s).*(?=\sseconds)','match'));
tline3=fgetl(fileID);

dataArray = textscan(fileID, formatspec , 'Delimiter', ',');
time=[char(dataArray{5}), repmat(' ',size(dataArray{1},1)), char(dataArray{6})];
time=datenum(time);
temp=dataArray{1,1}; %deg C
conductivity=dataArray{1,2}; %Siemens/m (S/m)
pres=dataArray{1,3};
sal=dataArray{1,4}; %salinity 

SB37.temperature.data=temp;
SB37.temperature.units='deg C';
SB37.conductivity.data=conductivity;
SB37.conductivity.units='S/m';
SB37.pressure.data=pres;
SB37.pressure.units='dB';
SB37.salinity.data=sal;
Sb37.salinity.units='g/kg'
SB37.time=time;

figure(); clf
subplot(411);
    plot(SB37.time,SB37.pressure.data, 'Color', 'green', 'LineWidth',1.5);
    ylabel('pres (dB)'); title('Raw Data'); set(gca, 'XTickLabel', []);
subplot(412);
    plot(SB37.time, SB37.temperature.data, 'Color', 'black','LineWidth',1.5);
    ylabel('temp (deg C)'); set(gca, 'XTickLabel', []);
subplot(413);
    plot(SB37.time, SB37.conductivity.data,'Color','blue', 'LineWidth',1.5);
    ylabel('cond (S/m)'); set(gca, 'XTickLabel', []);
subplot(414);
    plot(SB37.time, SB37.salinity.data,'Color', 'magenta' , 'LineWidth',1.5); 
    ylabel('sal (g/kg)');
    
   




