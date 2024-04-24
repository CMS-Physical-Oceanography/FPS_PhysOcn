% This file processes all raw instrument files to L0 MATFILES:
% Removes out of water portions
% Adds basic QC flags
% Runs on the BOEM OneDrive raw files folders
% Make a new file for each station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% 1. Add the local path to all the BOEM functions:
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/MacLawhorn/BOEM_Dep1/FunctionNToolbox/'))

% 2. Inputs:
% Deployment number (1 - 8):
depnum = 1;
% Station ID (should match folder names):
StationID='FPSC0';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/Documents/BOEM_OD/data/BOEM_deployment1';
% Exported file names for each instrument to process:
fname_c6 = 'FPSC_C6.csv';
fname_sbe = 'FPS-C_DEP1.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert:
start_time = '10/09/2023 16:00:00';
end_time = '10/30/2023 14:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath, filesep, StationID,filesep];
save_to_path = [deploypath,filesep,'L0',filesep];

% Instrument path names for the deployment. Change and add to as needed
fname_c6_full = [deploypath,'C6',filesep,fname_c6];
fname_sbe_full = [deploypath,'SBE37',filesep,fname_sbe];

%3. Run the processing codes:
%[C6L0,fname_c6_new] = C6_load(fname_c6_full,save_to_path,depnum,StationID,stime,etime);
% NOTE, this SBE had an odd time-base, need to change by 1 day below:
[SBE37L0,fname_sbe_new]=SBE37_load_nopressure(fname_sbe_full,save_to_path,depnum,StationID,stime,etime+2);

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on C6 (see RAW plots produced):
clearvars -except fname_c6_new
load(fname_c6_new)
[cdomf,qcFlag.CDOM] = optics_filter(CDOM,0.99);
%[chlf,qcFlag.chlorophyll_a] = optics_filter(C6.chlorophyll_a,0.99);
[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
CDOM = cdomf;turbidity = turbf;
notesqc = {'Removed spikes in CDOM and TURB.';
    sprintf('Run on %s',datestr(now))};
save(fname_c6_new,'SN','time','CDOM','chlorophyll_a','turbidity','depth','temperature','units','notes','notesqc','qcFlag'); 

%  Also have to shift the time of the SBE
clearvars -except fname_sbe_new stime etime
load(fname_sbe_new)
time2 = time-1;
aa = find(time2>=stime & time2 <=etime);
time = time2(aa);
temperature = temperature(aa);
conductivity = conductivity(aa);
pressure = pressure(aa);
salinity = salinity(aa);
notesqc = {'Had to shift the time back exactly 24 hours - careful with this SN.';
    sprintf('Run on %s',datestr(now))};
save(fname_sbe_new,'SN','time','temperature','conductivity','salinity','pressure','units','notes'); 
