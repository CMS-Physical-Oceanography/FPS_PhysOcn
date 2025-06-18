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
depnum = 5;
% Station ID (should match folder names):
StationID='FPSS1';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep5/';
% Exported file names for each instrument to process:
fname_rbrd1 = '207222_20241219_1505.rsk';
fname_rbrs = '200063_20241219_1512.rsk';
fname_sbe = '3571.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard and daylight are different:
start_time = '11/19/2024 16:00:00';
end_time = '12/09/2024 16:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath, filesep, StationID,filesep];
save_to_path = [deploypath,filesep,'L0',filesep];

%%
% Instrument path names for the deployment. Change and add to as needed
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_rbrs_full = [deploypath,fname_rbrs];
fname_sbe_full = [deploypath,fname_sbe];


%%
%3. Run the processing codes (better one at a time):
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[RBRs,fname_rbrs_new] = RBRsolot_load(fname_rbrs_full,save_to_path,depnum,StationID,stime,etime);
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime);

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc: