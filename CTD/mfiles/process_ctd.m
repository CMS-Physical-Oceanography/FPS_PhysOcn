% "rosette", "cape-fear", "both", "seahawk", "flow-through" or "all"?
ctdType = 'seahawk';
% ctdType = 'rosette';
% ctdType = 'cape-fear';
%
% data are stored in sub-dirs based on
% sample date: yyyymmdd
% dates = {'20230425'};%; %{'20230425'};
dates = {'20230619'};
% dates previously processed:
% ???????
%
% raw-data parent directory
[stat,rootDir] = system('git rev-parse --show-toplevel');
rootDir= rootDir(1:end-1);% drop the \n character
rawDir = [rootDir,filesep,'CTD/data/'];
% where to archive data?
arcDir = [rootDir,filesep,'CTD/mat_data/'];
%
cfDir = 'capeFear_flowThrough/';
shDir = 'seahawk_flowThrough/';    
rsDir = 'rosette/';
%
switch ctdType
  case 'cape-fear'
    datDir = [rawDir,cfDir];
    data = convert_capeFear_flowThrough_hex2mat(datDir,dates,arcDir);
  case 'rosette'
    datDir = [rawDir,rsDir];
    data = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);    
  case 'both'
    datDir = [rawDir,cfDir];
    cfData = convert_capeFear_flowThrough_hex2mat(datDir,dates,arcDir);
    datDir = [rawDir,rsDir];
    rsData = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);
  case 'seahawk'
    datDir = [rawDir,shDir];
    data = convert_seaHawk_flowThrough_to_mat(datDir,dates,arcDir);
end


