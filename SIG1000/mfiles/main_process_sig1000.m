clear all
close all
% stages of processing
% 1) define deployment number:
deploy  = 1
% 2) raw data input directory & filename convention:
rootDIR = ['/Users/derekgrimes/OneDriveUNCW/DATA/BOEM/FPSD1_Sig1000/mat_data/'];
fRoot   = ['S103071A010_FPS1_'];
% 3) output directory:
outRoot = ['/Users/derekgrimes/OneDriveUNCW/SIG1000/'];
% 4) output data file prefix:
filePrefix= sprintf('SIG_00103071_DEP%d_FPSC0_',deploy);
L0FRoot  = sprintf('%sL0_',filePrefix);
%
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
L1FRoot   = sprintf('%sL1',filePrefix);
%
% 5) time-periods when instrument was air (leave times empty to manually reselect them)
atmosphTime = [datenum(['09-Oct-2023 12:28:13'; '09-Oct-2023 15:12:12'])';...
               datenum(['30-Oct-2023 15:01:26'; '31-Oct-2023 04:35:57'])'];
% 6) deploy/recovery times
deployTime  = datenum('09-Oct-2023 16:00:00');
recoverTime = datenum('30-Oct-2023 14:00:00');
%
files = dir([rootDIR,fRoot,'*.mat']);
Nf    = length(files);
%
% height of transducer off of the bottom;
hab = 0.5;
HeadingOffset = -9;
%
% read in first/last data files and estimate atmospheric pressures
select_times_from_pressure
%
% load and pre-process data.
load_and_processes_sig1000_matrix_format
%
% make time-averages
time_average_and_rotate_sig1000
%
% estimate hourly wave stats
wave_statistics_sig1000