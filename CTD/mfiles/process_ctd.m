% "rosette", "cape-fear", "both", "seahawk", "flow-through" or "all"?
ctdType = 'seahawk';
% ctdType = 'rosette';
%
% data are stored in sub-dirs based on
% sample date: yyyymmdd
% dates = {'20230425'};%; %{'20230425'};
dates = {'20230619'};
% dates previously processed:
% ???????
%
% raw-data parent directory
rawDir = '/Users/derekgrimes/Projects/FPS/data/';
% where to archive data?
arcDir = '/Users/derekgrimes/Projects/FPS/mat_data/';
%
cfDir = 'capeFear_flowThrough/';
shDir = 'seahawk_flowThrough/';    
rsDir = 'rosette/';
%
switch ctdType
  case 'cape-fear'
    datDir = [rawDir,cfDir];
    SBE = convert_capeFear_flowThrough_hex2mat(datDir,dates,arcDir);
  case 'rosette'
    datDir = [rawDir,rsDir];
    SBE = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);    
  case 'both'
    datDir = [rawDir,cfDir];
    cfSBE = convert_capeFear_flowThrough_hex2mat(datDir,dates,arcDir);
    datDir = [rawDir,rsDir];
    rsSBE = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);
  case 'seahawk'
    datDir = [rawDir,shDir];
    SBE = convert_seaHawk_flowThrough_to_mat(datDir,dates,arcDir);
end
