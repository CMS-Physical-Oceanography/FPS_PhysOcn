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
depnum = 2;
% Station ID (should match folder names):
StationID='FPSS1';
% Base path to OneDrive or wherever the data files lie:
%basepath = '~/Documents/BOEM_OD/data/BOEM_deployment2';
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep2/';
% Exported file names for each instrument to process:
fname_rbrt = '214016_20240410_2055.rsk';
fname_rbrd = '207221_20240410_1944.rsk';
fname_sbe = 'SBE37_CTD_FPSS1_dep2.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard and daylight are different:
start_time = '02/15/2024 18:00:00';
end_time = '04/08/2024 17:00:00';

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
fname_rbrt_full = [deploypath,'RBR',filesep,fname_rbrt];
fname_rbrd_full = [deploypath,'RBR',filesep,fname_rbrd];
fname_sbe_full = [deploypath,'SBE37',filesep,fname_sbe];

%%
%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime); % 11 hour time difference??? Need to look at pressure
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd_full,save_to_path,depnum,StationID,stime,etime);
[RBRt,fname_rbrt_new] = RBRtri_load(fname_rbrt_full,save_to_path,depnum,StationID,stime,etime);

%% NEED TO DO ADCP!
% ADCP - first create a structure that includes a few essentials
ADCPparams.surfrange = [15:25]; % the approximate range of bins where we expect the surface to be
ADCPparams.bins2del = 2; % the number of bins below the surface that will be masked
ADCPparams.beamnum = 1; % the beam to use to estimate the surface from intensity
ADCPparams.qctouse = 'intn'; % the field to mask to use for QC
ADCPparams.mounting_dist = 0.5; % the approximate distance from the seabed.
% Currents
inpathc = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep2/FPSS1/RDIWH/WAVES/WAVES_000_000_TS2402151300_CUR.PD0';
[currents] = RDIWH_load_currents(inpathc,save_to_path,depnum,StationID,ADCPparams,stime,etime);
% Waves path to directory
inpathw = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep2/FPSS1/RDIWH/WAVES/';
[waves] = RDIWH_load_waves(inpathw,save_to_path,depnum,StationID,stime,etime);
filename=[save_to_path,'RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
save(filename,'currents','waves')

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% SBE needs to be time-shifted by 12 hours because we messed up the programming:
clearvars -except fname_sbe_new
load(fname_sbe_new)
time = time+0.5;
end_time = '04/08/2024 17:00:00';etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');
aa = find(time<=etime);
time = time(aa);temperature = temperature(aa);
conductivity = conductivity(aa);pressure = pressure(aa);
salinity = salinity(aa);
notesqc = {'Had to shift time by 12 hours - missing beginning';
    sprintf('Run on %s',datestr(now))};
save(fname_sbe_new,'SN','time','temperature','conductivity','salinity','pressure','units','notes','notesqc'); 


% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
clearvars -except fname_rbrt_new
load(fname_rbrt_new)
% All channels look bad after this time:
end_time = '03/22/2024 23:48:00';etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');
aa = find(time<=etime);
time = time(aa);chlorophyll_a = chlorophyll_a(aa);
FDOM = FDOM(aa);turbidity = turbidity(aa);
[fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.9);
%[chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.8);
FDOM = fdomf;turbidity = turbf;
notesqc = {'Removed spikes in FDOM and TURB. Data chopped after 03/22/24 23:48:00 because suspect';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 