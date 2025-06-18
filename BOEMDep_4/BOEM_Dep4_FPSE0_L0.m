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
StationID='FPSE0';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/OneDrive - UNC-Wilmington/Documents/Research/BOEM_OneDrive/General/data/BOEM_deployment4';
% Exported file names for each instrument to process:
fname_rbrs = '200065_20240922_1834.rsk';
fname_rbrd1 = '210872_20240922_1822.rsk';
fname_rbrt = '214016_20240922_2108.rsk';


fname_sbe = 'SN3571.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard and daylight are different:
start_time = '08/20/2024 18:00:00';
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

%%
% Instrument path names for the deployment. Change and add to as needed
fname_rbrs_full = [deploypath,fname_rbrs];
fname_rbrd1_full = [deploypath,fname_rbrd1];
fname_sbe_full = [deploypath,fname_sbe];
fname_rbrt_full = [deploypath,fname_rbrt];


%%
%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime);
[RBRs,fname_rbrs_new] = RBRsolot_load(fname_rbrs_full,save_to_path,depnum,StationID,stime,etime);
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd1_full,save_to_path,depnum,StationID,stime,etime);
[RBRt,fname_rbrt_new] = RBRtri_load(fname_rbrt_full,save_to_path,depnum,StationID,stime,etime);

%% ADCP
clearvars -except basepath save_to_path depnum StationID stime etime

% ADCP - first create a structure that includes a few essentials
ADCPparams.surfrange = [15:25]; % the approximate range of bins where we expect the surface to be
ADCPparams.bins2del = 2; % the number of bins below the surface that will be masked
ADCPparams.beamnum = 1; % the beam to use to estimate the surface from intensity
ADCPparams.qctouse = 'intn'; % the field to mask to use for QC
ADCPparams.mounting_dist = 0.5; % the approximate distance from the seabed.
% Currents
inpathc = strcat(basepath,'/',StationID,'/RAW/WAVES/WAVES_000_000_CUR.PD0');
[currents] = RDIWH_load_currents(inpathc,save_to_path,depnum,StationID,ADCPparams,stime,etime);
% Waves path to directory
inpathw = strcat(basepath,'/',StationID,'/RAW/WAVES/SPEC/');
[waves] = RDIWH_load_waves(inpathw,save_to_path,depnum,StationID,stime,etime);
filename=[save_to_path,'RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
save(filename,'currents','waves')

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
clearvars -except fname_rbrt_new
load(fname_rbrt_new)
[fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.99);
[chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99); %this time around its real bad:
aa = find(chlf>=1e4); 
chlf(aa) = NaN;qcFlag.chlorophyll_a(aa) = 0;
[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.9);
FDOM = fdomf;turbidity = turbf;
chlorophyll_a = chlf;
notesqc = {'Removed spikes in FDOM and TURB and Chlf. The end-time could be chopped';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 