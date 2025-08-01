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
StationID='FPSE1';
% Base path to OneDrive or wherever the data files lie:
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep6/';
% Exported file names for each instrument to process:
fname_rbrd1 = '210871_20250401_1823.rsk';
fname_rbrd2 = '210872_20250401_1838.rsk';
fname_rbrt = '214015_20250402_1731.rsk';
%fname_sbe = '3571.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs.
start_time = '02/04/2025 16:00:00';
end_time = '03/28/2025 15:00:00';

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
fname_rbrd2_full = [deploypath,fname_rbrd2];
fname_tri_full = [deploypath,fname_rbrt];

fname_sbe_full = [deploypath,fname_sbe];


%%
%3. Run the processing codes - likely have to be in local drive until we can figure out the RBR issue:
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd2_full,save_to_path,depnum,StationID,stime,etime);
[RBRt,fname_rbrt_new] = RBRtri_load(fname_tri_full,save_to_path,depnum,StationID,stime,etime);


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

%% ADCP - this is a REAL pain:
% ADCP - first create a structure that includes a few essentials
%cd(basepath)
ADCPparams.surfrange = [10:20]; % the approximate range of bins where we expect the surface to be
ADCPparams.bins2del = 2; % the number of bins below the surface that will be masked
ADCPparams.beamnum = 1; % the beam to use to estimate the surface from intensity
ADCPparams.qctouse = 'intn'; % the field to mask to use for QC
ADCPparams.mounting_dist = 0.5; % the approximate distance from the seabed.

[ADCP1] = rdradcp('WAVES/WAVES_000_000_TS2502041400_CUR.PD0',1,-1);
[ADCP1q] = edit_ADCP_currents_dep6(ADCP1,ADCPparams);
[ADCP2] = rdradcp('WAVESv2/WAVESv3_000_000_TS2502091930_CUR.PD0',1,-1);
[ADCP2q] = edit_ADCP_currents_dep6(ADCP2,ADCPparams);
[ADCP3] = rdradcp('WAVESv3/WAVESv3_000_000_TS2502201935_CUR.PD0',1,-1);
[ADCP3q] = edit_ADCP_currents_dep6(ADCP3,ADCPparams);
[ADCP4] = rdradcp('WAVESv4/WAVESv4_000_000_TS2502222205_CUR.PD0',1,-1);
[ADCP4q] = edit_ADCP_currents_dep6(ADCP4,ADCPparams);
[ADCP5] = rdradcp('WAVESv5/WAVESv5_000_000_TS2502251250_CUR.PD0',1,-1);
[ADCP5q] = edit_ADCP_currents_dep6(ADCP5,ADCPparams);
[ADCP6] = rdradcp('WAVESv6/WAVESv6_000_000_TS2502262015_CUR.PD0',1,-1);
[ADCP6q] = edit_ADCP_currents_dep6(ADCP6,ADCPparams);
[ADCP7] = rdradcp('WAVESv7/WAVESv7_000_000_TS2503011120_CUR.PD0',1,-1);
[ADCP7q] = edit_ADCP_currents_dep6(ADCP7,ADCPparams);
[ADCP8] = rdradcp('WAVESv8/WAVESv8_000_000_TS2503021845_CUR.PD0',1,-1);
[ADCP8q] = edit_ADCP_currents_dep6(ADCP8,ADCPparams);
[ADCP9] = rdradcp('WAVESv9/WAVESv9_000_000_TS2503040210_CUR.PD0',1,-1);
[ADCP9q] = edit_ADCP_currents_dep6(ADCP9,ADCPparams);
[ADCP10] = rdradcp('WAVESv10/WAVESv10_000_000_TS2503050930_CUR.PD0',1,-1);
[ADCP10q] = edit_ADCP_currents_dep6(ADCP10,ADCPparams);

% Run concatenate_ADCP_currents_dep6.m
filename=['L0/RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
% STILL NEED TO RUN WAVES!
save(filename,'currents')


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