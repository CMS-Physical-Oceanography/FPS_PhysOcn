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
depnum = 7;
% Station ID (should match folder names):
StationID='FPSC0';
% Base path to OneDrive or wherever the data files lie:
basepath = '~/OneDrive - UNC-Wilmington/Documents/Research/BOEM_OneDrive/General/data/BOEM_deployment7';

% Exported file names for each instrument to process:
fname_rbrd1 = '210870_20250604_1827.rsk';
fname_rbrd2 = '210872_20250618_1446.rsk';
fname_rbrc = '210866_20250604_1834.rsk';
fname_sbe = '1311.asc';
fname_c6 = 'C6.csv';

% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard (what it is in winter GMT-5) and daylight (what it is now GMT-4) are different:
start_time = '05/08/2025 16:00:00';
end_time = '06/03/2025 15:00:00';

% ------------------------
% ENDs TYPICAL USER INPUT besides the *_load.m commands (see 3. Below)
% ------------------------
% Make MATLAB datenum:
stime = datenum(start_time,'mm/dd/yyyy HH:MM:SS');
etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');

% Creates standard paths we use:
deploypath=[basepath,filesep, StationID,filesep,'RAW',filesep];
save_to_path = [basepath,filesep, StationID,filesep,'L0',filesep];
%%
% Instrument path names for the deployment. Change and add to as needed
fname_sbe_full = [deploypath,fname_sbe];
fname_c6_full = [deploypath,fname_c6];
fname_ctd_full = [deploypath,fname_sbe];

%%
%3. Run the processing codes - REMEMBER, these RBR codes can't seem to run on the remote drive:
[RBRctd, filename] = RBRctd_load(fname_rbrc,save_to_path,depnum,StationID,stime,etime);
[RBRd1,fname_rbrd_new] = RBRduet_load(fname_rbrd1,save_to_path,depnum,StationID,stime,etime);
[RBRd2,fname_rbrd_new] = RBRduet_load(fname_rbrd2,save_to_path,depnum,StationID,stime,etime);

[SBE37L0,fname_sbe_new]=SBE37_load_nopressure(fname_sbe_full,save_to_path,depnum,StationID,stime,etime); 
[C6L0,fname_c6_new] = C6_load_ALL(fname_c6_full,save_to_path,depnum,StationID,stime,etime);

%% ADCP
clearvars -except basepath save_to_path depnum StationID stime etime

% ADCP - first create a structure that includes a few essentials
ADCPparams.surfrange = [15:25]; % the approximate range of bins where we expect the surface to be
ADCPparams.bins2del = 2; % the number of bins below the surface that will be masked
ADCPparams.beamnum = 1; % the beam to use to estimate the surface from intensity
ADCPparams.qctouse = 'intn'; % the field to mask to use for QC
ADCPparams.mounting_dist = 0.5; % the approximate distance from the seabed.
% Currents
inpathc = strcat(basepath,'/',StationID,'/RAW/ADCP/WAVES_000_000_TS2505081100_CUR.PD0');
[currents] = RDIWH_load_currents(inpathc,save_to_path,depnum,StationID,ADCPparams,stime,etime);
% Waves path to directory
inpathw = strcat(basepath,'/',StationID,'/RAW/ADCP/SPEC/');
[waves] = RDIWH_load_waves(inpathw,save_to_path,depnum,StationID,stime,etime);
filename=[save_to_path,'RDI_',sprintf('%08d',currents.SN),'_','DEP',num2str(depnum),'_',StationID,'_L0.mat'];
save(filename,'currents','waves')

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% In this case, we do some spike filtering on RBRtri (see RAW plots produced):
fname_sbe_new = [save_to_path,'SBE_00001311_DEP7_FPSC0_L0.mat'];
fname_c6_new = [save_to_path,'0C6_23600142_DEP7_FPSC0_L0.mat'];
fname_ctd_new = [save_to_path,'RBR_00210866_DEP7_FPSC0_L0.mat'];

clearvars -except fname_c6_new fname_sbe_new fname_ctd_new
load(fname_sbe_new)
% [cdomf,qcFlag.FDOM] = optics_filter(CDOM,0.99);
% [chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.99);
 [turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
% FDOM = fdomf;turbidity = turbf;
% chlorophyll_a = chlf;
% notesqc = {'Removed spikes in FDOM and TURB and Chlf.';
%     sprintf('Run on %s',datestr(now))};
% save(fname_rbrt_new,'SN','time','FDOM','chlorophyll_a','turbidity','units','notes','notesqc','qcFlag'); 