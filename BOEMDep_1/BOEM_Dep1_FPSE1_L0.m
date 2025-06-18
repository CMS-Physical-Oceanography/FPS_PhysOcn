% This file processes all raw instrument files to L0 MATFILES:
% Removes out of water portions
% Adds basic QC flags
% Runs on the BOEM OneDrive raw files folders
% Make a new file for each station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% 1. Add the local path to all the BOEM functions:
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/MacLawhorn/BOEM_Dep1/FunctionNToolbox/'))
addpath('../FunctionNToolbox')

% 2. Inputs:
% Deployment number (1 - 8):
depnum = 1;
% Station ID (should match folder names):
StationID='FPSE1';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/OneDrive - UNC-Wilmington/Documents/Research/BOEM_OneDrive/General/data/BOEM_deployment1';
% Exported file names for each instrument to process:
fname_rbrt = '214016_20231031_2014.rsk';
fname_rbrs = '200061_20231031_1951.rsk';
fname_sbe = 'FPS-E1_DEP1.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert:
start_time = '10/09/2023 18:00:00';
end_time = '10/30/2023 13:00:00';

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
fname_rbrt_full = [deploypath,'RBRtri',filesep,fname_rbrt];
fname_rbrs_full = [deploypath,'RBRsolot',filesep,fname_rbrs];
fname_sbe_full = [deploypath,'SBE37',filesep,fname_sbe];

%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,0,stime,etime);
[RBRs,fname_rbrs_new] = RBRsolot_load('FPSE1/200061_20231031_1951.rsk',save_to_path,depnum,StationID,stime,etime);
[RBRt,fname_rbrt_new] = RBRtri_load('FPSE1/214016_20231031_2014.rsk',save_to_path,depnum,StationID,stime,etime);

% ADCP - remember, bad pressure and bad waves
ADCPparams.surfrange = [15:25];
ADCPparams.bins2del = 2;
ADCPparams.beamnum = 1;
ADCPparams.qctouse = 'intn';
ADCPparams.mounting_dist = 0.5; 
inpath = strcat(basepath,'/',StationID,'/RDI_WH/WAVES/WAVES_000_000_TS2310091200_CUR.PD0')
[currents] = RDIWH_load_currents(inpath,save_to_path,depnum,StationID,ADCPparams,stime,etime);
filename=[save_to_path,'RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
save(filename,'currents')

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% The SBE has has some salinity spiking. Delete them:
clearvars -except fname_sbe_new stime etime
load(fname_sbe_new)
[salinity,qcFlag.salinity] = optics_filter(salinity,0.995);
[conductivity,qcFlag.conductivity] = optics_filter(conductivity,0.99);
notesqc = {'Conductivity and salinity spikes removed';
    sprintf('Run on %s',datestr(now))};
save(fname_sbe_new,'SN','time','temperature','conductivity','salinity','pressure','units','notes','notesqc','qcFlag'); 


% In this case, we do some spike filtering on RBR (see RAW plots produced):
clearvars -except fname_rbrt_new
load(fname_rbrt_new)
[fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.99);
[chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
FDOM = fdomf;chlorophyll_a = chlf;turbidity = turbf;
notesqc = {'Removed spikes in FDOM, ChL, and TURB.';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 