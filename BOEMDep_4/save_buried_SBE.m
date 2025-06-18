fname = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep4/FPSC0/BOEM4_26944_C0.asc';
% Process the new SBE - data file has some issues:
fileID=fopen(fname,'r');

    formatspec='%f %f %f %f %s %s';
   
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

% Read each line, skip th <execute> lines:
count = 1;
while ~feof(fileID)
    line = fgetl(fileID);

     is=findstr(line,',');
     if length(is) == 0 
     else
         DATA(count,1) = str2num(line(1:is(1)-1));
         DATA(count,2) = str2num(line(is(1)+1:is(2)-1));
         DATA(count,3) = str2num(line(is(2)+1:is(3)-1));
         DATA(count,4) = str2num(line(is(3)+1:is(4)-1));
         time = char(strcat(line(is(4)+1:is(5)-1),' ',line(is(5)+1:end)));
         TIME(count) = datenum(time);
         count = count + 1;
     end
end


SBE37.temperature=DATA(:,1);
SBE37.conductivity=DATA(:,2);
SBE37.pressure=DATA(:,3);
SBE37.salinity=DATA(:,4);
SBE37.time=TIME';
SBE37.SN= 26944;
SBE37.units={'temp = deg C'; 'cond = S/m'; 'pres = dB'; 'sal = PSU'};
SBE37.notes = {'Time base is GMT'; 'This was a buried sensor. The bad periods not removed. Note, time is non-uniform'; 
    sprintf('Created on %s',datestr(now))};

savename = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep4/FPSC0/L0/SBE_00026944_DEP4_FPSC0_L0.mat';
save(savename,'SBE37')