%This file processes all raw instrument files to L0 MATFILES:
% 1. Add the path to all the BOEM functions:
addpath(genpath('/Users/suandas/Documents/Research_Local/BOEM/BOEM_Dep1/MacLawhorn/BOEM_Dep1/FunctionNToolbox/'))

%%. 2. Repeat the following for each Station ID:

%% FPSC:
clear; clc;
depnum = 1; %Change as needed
basepath = '/Users/suandas/Documents/Research_Local/BOEM';%Change as needed
StationID='FPSC0';%Change as needed
deployname= strcat('BOEM_Dep',num2str(depnum));
deploypath=[basepath, filesep, deployname, filesep, StationID, filesep];
save_to_path = [StationID,filesep,'L0processing',filesep];

% Instrument path names for the deployment. Change and add to as needed
fnamec6 = [deploypath,'C6',filesep,'FPSC_C6.csv'];
fnamesbe = [deploypath,'SBE37',filesep,'FPS-C_DEP1.asc'];
% END user input
%% Run the codes:
[C6L0] = C6_load(fnamec6,save_to_path,depnum,StationID);
[SBE37L0]=SBE37_load_nopressure(fnamesbe,save_to_path,depnum,StationID);

%% FPSS1:
clear; clc;
depnum = 1;
basepath = '/Users/suandas/Documents/Research_Local/BOEM';
StationID='FPSS1';
deployname= strcat('BOEM_Dep',num2str(depnum));
deploypath=[basepath, filesep, deployname, filesep, StationID, filesep];
% Save path 
save_to_path = [StationID,filesep,'L0processing',filesep];

% Instrument path names for the deployment
fnamesbe = [deploypath,'SBE37',filesep,'FPS-S1_DEP1.asc'];
fnamerbrd = [deploypath,'207221_20231031_1958.rsk'];
fnamerbrt = [deploypath,'214015_20231031_2010.rsk'];
fnamewh = deploypath;

% END user input
%% Run the codes:
[SBE37L0]=SBE37_load(fnamesbe,save_to_path,depnum,StationID,1);
[RBRduetL0] = RBRduet_load(fnamerbrd,save_to_path,depnum,StationID);
[RBRtriL0] = RBRtri_load(fnamerbrt,save_to_path,depnum,StationID);
[RDIWHL0]=RDIWH_load(fnamewh,save_to_path,depnum,StationID)

%% FPSE1:
clear; clc;
depnum = 1;
basepath = '/Users/suandas/Documents/Research_Local/BOEM';
StationID='FPSE1';
deployname= strcat('BOEM_Dep',num2str(depnum));
deploypath=[basepath, filesep, deployname, filesep, StationID, filesep];
% Save path 
save_to_path = [StationID,filesep,'L0processing',filesep];

% Instrument path names for the deployment
fnamesbe = [deploypath,'SBE37',filesep,'FPS-E1_DEP1.asc'];
fnamerbrs = [deploypath,'200061_20231031_1951.rsk'];
fnamerbrt = [deploypath,'214016_20231031_2014.rsk'];
fnamewh = deploypath;

% END user input
%% Run the codes:
[SBE37L0]=SBE37_load(fnamesbe,save_to_path,depnum,StationID,0);
[RBRsolotL0] = RBRsolot_load(fnamerbrs,save_to_path,depnum,StationID);
[RBRtriL0] = RBRtri_load(fnamerbrt,save_to_path,depnum,StationID);
[RDIWHL0]=RDIWH_load(fnamewh,save_to_path,depnum,StationID) % LAST ONE!
