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
depnum = 3;
% Station ID (should match folder names):
StationID='FPSC0';
% Base path to OneDrive or wherever the data files lie:
%basepath = '~/Documents/BOEM_OD/data/BOEM_deployment3';
basepath = '/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep3/';
% Exported file names for each instrument to process:
fname_rbrd = '210870_20240606_1554.rsk';
fname_rbrs1 = '200065_20240606_1613.rsk';
fname_rbrs2 = '200064_20240606_1524.rsk';
fname_sbe = 'BOEM3_1311_C0.asc';
fname_c6 = 'C6BOEM3.csv';


% Start time and end time, rounded to nearsest hour that instrument is stable.
% In GMT - look at the deploy/recover logs and convert CAREFUL standard and daylight are different:
start_time = '05/13/2024 17:00:00';
end_time = '06/03/2024 15:00:00';

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
fname_rbrd_full = [deploypath,fname_rbrd];
fname_rbrs1_full = [deploypath,fname_rbrs1];
fname_rbrs2_full = [deploypath,fname_rbrs2];
fname_sbe_full = [deploypath,fname_sbe];
fname_c6_full = [deploypath,fname_c6];

%%
%3. Run the processing codes:
[SBE37L0,fname_sbe_new]=SBE37_load_nopressure(fname_sbe_full,save_to_path,depnum,StationID,stime,etime); % 11 hour time difference??? Need to look at pressure
[RBRs1,fname_rbrs_new] = RBRsolot_load(fname_rbrs1_full,save_to_path,depnum,StationID,stime,etime);
[RBRs2,fname_rbrs_new] = RBRsolot_load(fname_rbrs2_full,save_to_path,depnum,StationID,stime,etime);
[RBRd,fname_rbrd_new] = RBRduet_load(fname_rbrd_full,save_to_path,depnum,StationID,stime,etime);
[C6L0,fname_c6_new] = C6_load(fname_c6_full,save_to_path,depnum,StationID,stime,etime);

%% ADCP

%% 4. OPTIONALS - this will depend on the instrument, deployment, etc:
% % In this case, we do some spike filtering on C6 (see RAW plots produced):
% clearvars -except fname_c6_new
% load(fname_c6_new)
% [cdomf,qcFlag.CDOM] = optics_filter(CDOM,0.99);
% [chlf,qcFlag.chlorophyll_a] = optics_filter(chlorophyll_a,0.9);
% [turbf,qcFlag.turbidity] = optics_filter(turbidity,0.99);
% %ODD turbidity for the first 900 samples: 
% turbf(1:900) = NaN;
% [turbf2,qcFlag.turbidity] = optics_filter(turbf,0.9);
% 
% %Rename
% CDOM = cdomf;turbidity = turbf2;chlorophyll_a = chlf;
% notesqc = {'Removed spikes in CDOM, CHL and TURB. TUB first 900 data points also removed (suspect)';
%     sprintf('Run on %s',datestr(now))};
% save(fname_c6_new,'SN','time','CDOM','chlorophyll_a','turbidity','depth','temperature','units','notes','notesqc','qcFlag'); 
% 
% % SBE needs to be time-shifted by 12 hours because we messed up the programming:
% clearvars -except fname_sbe_new
% load(fname_sbe_new)
% time = time+0.5;
% end_time = '04/08/2024 16:00:00';etime = datenum(end_time,'mm/dd/yyyy HH:MM:SS');
% aa = find(time<=etime);
% time = time(aa);temperature = temperature(aa);
% conductivity = conductivity(aa);pressure = pressure(aa);
% salinity = salinity(aa);
% notesqc = {'Had to shift time by 12 hours - missing beginning';
%     sprintf('Run on %s',datestr(now))};
% save(fname_sbe_new,'SN','time','temperature','conductivity','salinity','pressure','units','notes','notesqc'); 
