%BOEM Deployment Test 1 Function Test File
clear; clc;

deployname= 'BOEMTest1_062923_063023';
crateNum='Crate1';
deploypath=['C:\Users\amcla\OneDrive - UNC-Wilmington\Documents\FPS\DeployData', filesep, deployname, filesep, crateNum];
addpath(genpath('C:\Users\amcla\OneDrive - UNC-Wilmington\Documents\Functions'));

[RDIWHL0]=RDIWH_load(deploypath)

[RBRtriL0]=RBRtri_load(deploypath)

[SBE37L0]=SBE37_load(deploypath)

[RBRsoloTL0]=RBRsoloT_load(deploypath)

[ADCP, RBRtri ,SBE37, RBRsoloT]=L1_BOEMprocessing(deploypath, 13, 26);

[L1Plots]plot_L1_QC(deploypath);

