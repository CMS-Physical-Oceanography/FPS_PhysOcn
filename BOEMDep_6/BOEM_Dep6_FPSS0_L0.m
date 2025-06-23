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
StationID='FPSS0';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep6/';
% Exported file names for each instrument to process:
fname_rbrd1 = '210869_20250421_2020.rsk';
fname_rbrs = '200061_20250421_2015.rsk';
fname_sbe = '3661.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs.
start_time = '02/04/2025 18:00:00';
end_time = '04/17/2025 14:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath, filesep, StationID,filesep];
save_to_path = [basepath, filesep, StationID,filesep,'L0',filesep];

%%
% Instrument path names for the deployment. Change and add to as needed
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_solo_full = [deploypath,fname_rbrs];

fname_sbe_full = [deploypath,fname_sbe];


%%
%3. Run the processing codes - likely have to be in local drive until we can figure out the RBR issue:
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[RBRs,fname_rbrs_new] = RBRsolot_load(fname_solo_full,save_to_path,depnum,StationID,stime,etime);
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime); 

%% The SBE used the new CTD. Don't have this as part of the standard
% processing yet. Convert to CNV through SBE-Data-Processing:
[data,names]=cnv2mat_BOEM('37-SMP_03726944_2025_04_02.cnv');
% Create the fields we want to save:
time = datenum(2025,1,data(:,1));
SN = '26944';
conductivity = data(:,2);
temperature = data(:,3);
pressure = data(:,5);
salinity = data(:,4);
units={'temp = deg C'; 'cond = mS/cm'; 'pres = dB'; 'sal = PSU'};
notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
% Trim the beginning and ends of the files:
aa = find(time>=stime & time <=etime);
time = time(aa);
temperature = temperature(aa);
conductivity = conductivity(aa);
pressure = pressure(aa);
salinity = salinity(aa);

save('L0/SBE_00026944_DEP6_FPSE1_L0.mat','SN','time','temperature','conductivity','salinity','pressure','units','notes'); 

%%
% %% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% % In this case, we do some spike filtering on RBRtri (see RAW plots produced):
% clearvars -except fname_rbrt_new
% load(fname_rbrt_new)
% [fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.99);
% [chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
% [turbf,qcFlag.turbidity] = optics_filter(turbidity,0.9);
% FDOM = fdomf;turbidity = turbf;
% chlorophyll_a = chlf;
% notesqc = {'Removed spikes in FDOM and TURB and Chlf.';
%     sprintf('Run on %s',datestr(now))};
% save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 