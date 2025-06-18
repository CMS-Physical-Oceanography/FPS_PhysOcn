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
StationID='FPSE0';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep6/';
% Exported file names for each instrument to process:
fname_rbrd1 = '210870_20250401_1849.rsk';
fname_sbe = '3571.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs.
start_time = '02/04/2025 17:00:00';
end_time = '03/28/2025 14:00:00';

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
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_sbe_full = [deploypath,fname_sbe];


%%
%3. Run the processing codes:
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime);
