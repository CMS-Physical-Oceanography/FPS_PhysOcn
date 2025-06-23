% This file processes all raw instrument files to L0 MATFILES:
% Removes out of water portions
% Adds basic QC flags
% Runs on the BOEM OneDrive raw files folders
% Make a new file for each station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% 1. Add the local path to all the BOEM functions:
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/MacLawhorn/BOEM_Dep1/FunctionNToolbox/'))
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/'))
addpath('../FunctionNToolbox')

% 2. Inputs:
% Deployment number (1 - 8):
depnum = 4;
% Station ID (should match folder names):
StationID='FPSE1';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep4/';
% Exported file names for each instrument to process:
%fname_rbrs = '200063_20240606_1620.rsk'; %Missing??
fname_rbrd1 = '210871_20240923_1405.rsk';
fname_c6 = 'BOEM4_C6.csv';

fname_sbe = 'SN3570.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard and daylight are different:
start_time = '08/20/2024 19:00:00';
end_time = '09/19/2024 14:00:00';

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
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_sbe_full = [deploypath,fname_sbe];
fname_c6_full = [deploypath,fname_c6];


%%
%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load_nopressure(fname_sbe_full,save_to_path,depnum,StationID,stime,etime); % 11 hour time difference??? Need to look at pressure
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[C6L0,fname_c6_new] = C6_load_ALL(fname_c6_full,save_to_path,depnum,StationID,stime,etime);

%% ADCP
%clearvars -except save_to_path depnum StationID stime etime


%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
clearvars -except fname_c6_new
load(fname_c6_new)
[cdomf,qcFlag.FDOM] = optics_filter(CDOM,0.99);
[chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
%[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.9);
CDOM = cdomf;
chlorophyll_a = chlf;
notesqc = {'Removed spikes in CDOM and Chlf. Maybe the end-time could be chopped?';
    sprintf('Run on %s',datestr(now))};
save(fname_c6_new,'SN','time','CDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 