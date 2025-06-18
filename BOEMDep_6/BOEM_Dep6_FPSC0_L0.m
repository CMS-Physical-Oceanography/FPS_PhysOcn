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

% 2. Inputs:
% Deployment number (1 - 8):
depnum = 6;
% Station ID (should match folder names):
StationID='FPSC0';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep6/';
% Exported file names for each instrument to process:
fname_rbrd1 = '207221_20250401_1904.rsk';
fname_rbrs = '200065_20250401_1835.rsk';
fname_sbe = '3570.asc';
fname_c6 = 'C6.csv';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard (what it is in winter GMT-5) and daylight (what it is now GMT-4) are different:
start_time = '02/04/2025 18:00:00';
end_time = '03/28/2025 13:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath,filesep, StationID,filesep,'RAW',filesep];
save_to_path = ['/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep6/FPSC0/','L0',filesep];

%%
% Instrument path names for the deployment. Change and add to as needed
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_rbrs_full = [deploypath,fname_rbrs];
fname_sbe_full = [deploypath,fname_sbe];
fname_c6_full = [deploypath,fname_c6];


%%
%3. Run the processing codes:
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[RBRs,fname_rbrs_new] = RBRsolot_load(fname_rbrs_full,save_to_path,depnum,StationID,stime,etime);
[SBE37L0,fname_sbe_new]=SBE37_load_nopressure(fname_sbe_full,save_to_path,depnum,StationID,stime,etime); 
[C6L0,fname_c6_new] = C6_load(fname_c6_full,save_to_path,depnum,StationID,stime,etime);

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
clearvars -except fname_c6_new
% load(fname_c6_new)
% [fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.99);
% [chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
% [turbf,qcFlag.turbidity] = optics_filter(turbidity,0.9);
% FDOM = fdomf;turbidity = turbf;
% chlorophyll_a = chlf;
% notesqc = {'Removed spikes in FDOM and TURB and Chlf.';
%     sprintf('Run on %s',datestr(now))};
% save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 