% This file processes all raw instrument files to L0 MATFILES:
% Removes out of water portions
% Adds basic QC flags
% Runs on the BOEM OneDrive raw files folders
% Make a new file for each station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% 1. Add the local path to all the BOEM functions:
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/MacLawhorn/BOEM_Dep1/FunctionNToolbox/'))

% 2. Inputs:
% Deployment number (1 - 8):
depnum = 1;
% Station ID (should match folder names):
StationID='FPSS1';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/Documents/BOEM_OD/data/BOEM_deployment1';
% Exported file names for each instrument to process:
fname_rbrt = '214015_20231031_2010.rsk';
fname_rbrd = '207221_20231031_1958.rsk';
fname_sbe = 'FPS-S1_DEP1.asc';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert:
start_time = '10/10/2023 16:00:00';
end_time = '10/30/2023 15:00:00';

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
fname_rbrd_full = [deploypath,'RBRduet',filesep,fname_rbrd];
fname_sbe_full = [deploypath,'SBE37',filesep,fname_sbe];

%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load(fname_sbe_full,save_to_path,depnum,StationID,1,stime,etime);
[RBRd,fname_rbrd_new] = RBRduet_load('FPSS1/207221_20231031_1958.rsk',save_to_path,depnum,StationID,stime,etime);
[RBRt,fname_rbrt_new] = RBRtri_load('FPSS1/214015_20231031_2010.rsk',save_to_path,depnum,StationID,stime,etime);

%% NEED TO DO ADCP!
%load('FPSS1/L0processing/RDI_00011548_DEP1_FPSS1_L0.mat')

% ADCP - first create a structure that includes a few essentials
ADCPparams.surfrange = [15:25]; % the approximate range of bins where we expect the surface to be
ADCPparams.bins2del = 2; % the number of bins below the surface that will be masked
ADCPparams.beamnum = 1; % the beam to use to estimate the surface from intensity
ADCPparams.qctouse = 'intn'; % the field to mask to use for QC
ADCPparams.mounting_dist = 0.5; % the approximate distance from the seabed.
% Currents
inpathc = 'FPSS1/RDI_WH/WAVES/WAVES_000_000_TS2310101200_CUR.PD0';
[currents] = RDIWH_load_currents(inpathc,save_to_path,depnum,StationID,ADCPparams,stime,etime);
% Waves path to directory
inpathw = 'FPSS1/RDI_WH/WAVES/SPEC/';
[waves] = RDIWH_load_waves(inpathw,save_to_path,depnum,StationID,stime,etime);
filename=[save_to_path,'RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
save(filename,'currents','waves')

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% The RBR Duet has some odd behavior at the end including the time. Delete them:
clearvars -except fname_rbrd_new stime etime
load(fname_rbrd_new)
aa = find(time>=datenum(2023,10,29,8,0,0));
time(aa) = [];
temperature(aa) = [];
pressure(aa) = [];
notesqc = {'The last day has been deleted due to odd sensor behavior';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrd_new,'SN','time','temperature','pressure','units','notes','notesqc'); 


% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
clearvars -except fname_rbrt_new
load(fname_rbrt_new)
[fdomf,qcFlag.FDOM] = optics_filter(FDOM,0.99);
%[chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
[turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
FDOM = fdomf;turbidity = turbf;
notesqc = {'Removed spikes in FDOM and TURB.';
    sprintf('Run on %s',datestr(now))};
save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 