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
depnum = 7;
% Station ID (should match folder names):
StationID='FPSE1';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/OneDrive - UNC-Wilmington/Documents/Research/BOEM_OneDrive/General/data/BOEM_deployment7';

% Exported file names for each instrument to process:
fname_rbrd1 = '207221_20250626_2202.rsk';
fname_rbrd2 = '207222_20250626_2143.rsk';
fname_rbrc = '210865_20250626_2152.rsk';
fname_rbrt = '214015_20250626_2200.rsk';
fname_sbe = '3571.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard (what it is in winter GMT-5) and daylight (what it is now GMT-4) are different:
start_time = '05/08/2025 14:00:00';
end_time = '06/26/2025 14:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath, filesep, StationID,filesep,'RAW',filesep];
save_to_path = [basepath, filesep, StationID,filesep,'L0',filesep];

%%
% Instrument path names for the deployment. Change and add to as needed
fname_sbe_full = [deploypath,fname_sbe];


%%
%3. Run the processing codes - likely have to be in local drive until we can figure out the RBR issue:
[RBRctd, filename] = RBRctd_load(fname_rbrc,save_to_path,depnum,StationID,stime,etime);
[RBRd1,fname_rbrd_new] = RBRduet_load(fname_rbrd1,save_to_path,depnum,StationID,stime,etime);
[RBRd2,fname_rbrd_new] = RBRduet_load(fname_rbrd2,save_to_path,depnum,StationID,stime,etime);

[RBRt,fname_rbrt_new] = RBRtri_load(fname_rbrt,save_to_path,depnum,StationID,stime,etime);
[SBE37L0,fname_sbe_new]=SBE37_load('3571.asc',save_to_path,depnum,StationID,1,stime,etime); 



%%
% %% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
fname_sbe_new = [save_to_path,'SBE_00003571_DEP7_FPSE1_L0.mat'];
fname_rbrt_new = [save_to_path,'RBR_00214015_DEP7_FPSE1_L0.mat'];
fname_ctd_new = [save_to_path,'RBR_00210865_DEP7_FPSE1_L0.mat'];

clearvars -except fname_rbrt_new fname_sbe_new fname_ctd_new
%load(fname_sbe_new)
 load(fname_rbrt_new)
 [fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.9);
 [chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
 [turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
FDOM = fdomf;turbidity = turbf;
chlorophyll_a = chlf;
notesqc = {'Removed spikes in FDOM and TURB and Chlf.';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 

%Unfortunately, the salinity is bad:
load(fname_ctd_new)
 [salinityf,qcFlag.salinity] = optics_filter(salinity,0.99);
salinityf(646159:end) = NaN;
qcFlag.salinity(646159:end) = 1;
salinity = salinityf;
notesqc = {'Salinity is incorrect -  most of it is removed.';
    sprintf('Run on %s',datestr(now))};
save(fname_ctd_new,'SN','time','temperature','pressure','conductivity','salinity','units','notes','notesqc','qcFlag'); 
