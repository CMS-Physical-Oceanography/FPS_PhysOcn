function [data,names]=cnv2mat_BOEM(cnv_file);
% cnv2mat_BOEM Reads the SeaBird ASCII .CNV file format
%
%  Usage:   [data,names]=cnv2matAS(cnv_file);
%
%     Input:  cnv_file = name of .CNV file  (e.g. 'cast002.cnv')
%
%     Output:
%            data = matrix containing all the columns of data in the .CNV file
%           names = string matrix containing the names and units of the columns
%         sensors = string matrix containing the names of the sensors (commented out)
%
%  NOTE: Doesn't get lat lon or much else from the header data. Just readings

% Open the .cnv file as read-only text
%
  fid=fopen(cnv_file,'rt');
% 
% Read the header.
% Start reading header lines of .CNV file,
% Stop at line that starts with '*END*'
%
% Pull out NMEA lat & lon along the way and look
% at the '# name' fields to see how many variables we have.
%
str='*START*';
while (~strncmp(str,'*END*',5));
     str=fgetl(fid);
%-----------------------------------
%
     if (strncmp(str,'# name',6))  
        var=sscanf(str(7:10),'%d',1);
        var=var+1;  % .CNV file counts from 0, Matlab counts from 1
 %      stuff variable names into cell array
      names{var}=str;
%------------------------------
%  pick up bad flag value
     elseif (strncmp(str,'# bad_flag',10))  
        isub=13:length(str);
        bad_flag=sscanf(str(isub),'%g',1);
     end
end
%==============================================
%
%  Done reading header.  Now read the data!
%
nvars=var;  %number of variables

% Read the data into one big matrix
%
data=fscanf(fid,'%f',[nvars inf]);

fclose(fid);

%
% Flag bad values with nan
%
ind=find(data==bad_flag);
data(ind)=data(ind)*nan;

%
% Flip data around so that each variable is a column
data=data.';

% Convert cell arrays of names to character matrices
names=char(names);
% Can be commented out if don't need sensor info:
%sensors=char(sensors);

return
