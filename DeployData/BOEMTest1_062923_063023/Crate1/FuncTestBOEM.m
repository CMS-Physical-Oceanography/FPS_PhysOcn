%BOEM Deployment Test 1 Function Test File
clear; clc;
addpath(genpath('C:\Users\amcla\OneDrive - UNC-Wilmington\Documents\FPS_BOEMProj\FunctionNToolbox'))

deployname= 'BOEMTest1_062923_063023';
crateNum='Crate1';
deploypath=['C:\Users\amcla\OneDrive - UNC-Wilmington\Documents\FPS_BOEMProj\DeployData', filesep, deployname, filesep, crateNum];

[RDIWHL0]=RDIWH_load(deploypath);

[RBRtriL0]=RBRtri_load(deploypath);

[SBE37L0]=SBE37_load(deploypath);

[RBRsoloTL0]=RBRsoloT_load(deploypath);

%Input variables to the L1BOEMprocessing function below:
deploydepth = 13; %depth at which the lander is deployed

surfacebin = 26; %It defines the sea surface bin # where the echo intensities 
% are > 90 associated with the water depth and bin size. The bin number can be determined
% by looking at the raw output from the RDIWH_load.m function in ADCP.cur.intens variable or 
% divide the mean depth during the deployment by the bin size to get a bin
% number at the surface approximately. The numeric bin # is important to set 
%a parameter for removing collected data above the sea surface

%The start and end time of each BOEM deployment should be stored in a field
%note book from the day of the deployment and retrieval date. Reference
%the notes for input.
startTime = datetime('06/29/2023 02:00:00', 'InputFormat', 'MM/dd/yyyy HH:mm:ss');
endTime = datetime('06/30/2023 09:55:00', 'InputFormat', 'MM/dd/yyyy HH:mm:ss');

[ADCP, RBRtri ,SBE37, RBRsoloT]=L1_BOEMprocessing(deploypath, deploydepth, surfacebin, startTime, endTime);

[L1Plots]=plot_L1_QC(deploypath);

