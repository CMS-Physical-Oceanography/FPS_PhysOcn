clear all
close all
% stages of processing
% 1) define deployment number:
deploy  = 3;
% need to shift time to convert EST to UTC
time_shift = 0/24;
% 2) raw data input directory & filename convention:
rootDIR = sprintf('/Users/derekgrimes/OneDriveUNCW/DATA/BOEM/FPSD%d_SIG1000/mat_data/',deploy);
fRoot   = sprintf('S103071A013_FPS%d_',deploy);
% 3) output directory:
outRoot = sprintf('/Users/derekgrimes/OneDriveUNCW/Documents-UNCW-BOEM-FryingPanShoals/General/data/BOEM_deployment%d/FPSC0/',deploy);
% 4) output data file prefix:
outDir  = [outRoot,filesep,'SIG1000',filesep];
filePrefix= sprintf('SIG_00103071_DEP%d_FPSC0_',deploy);
%
% 4a) current/echo average interval (seconds)
dtAvg     = 300;
L0dir     = [outRoot, filesep, 'L0',filesep];
L0FRoot   = sprintf('%sL0',filePrefix);
L1dir     = [outRoot, filesep, 'L1',filesep];
L1FRoot   = sprintf('%sL1',filePrefix);
%
% 5) time-periods when instrument was air (leave times empty to manually reselect them)
atmosphTime = [datenum('13-May-2024 11:30:54'), datenum('13-May-2024 15:46:17');...
               datenum('03-Jun-2024 16:39:11'), datenum('04-Jun-2024 02:41:20')];%[datenum('15-Feb-2024 13:12:05'), datenum('15-Feb-2024 14:08:32')];
% 6) deploy/recovery times
deployTime  = [];%[datenum('15-Feb-2024 15:00:00')]; %datenum('09-Oct-2023 16:00:00');
recoverTime = [];%[datenum('17-Mar-2024 15:00:00')]; %datenum('30-Oct-2023 14:00:00');
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
% shift time limits
atmosphTime = atmosphTime + time_shift;
deployTime  = deployTime + time_shift;
recoverTime = recoverTime + time_shift;
%
% load and pre-process data.
% load_and_processes_sig1000_matrix_format
load_and_process_sig1000_RDI_matrix_format
%
% make time-averages
time_average_and_rotate_sig1000_RDI_matrix_format
%
% estimate hourly wave stats
estimate_wave_bulk_stats_SIG1000_RDI_matrix_format
%
%
fig = figure;
imagesc(datetime(waves.Time','convertFrom','datenum'),waves.frequency,log10(waves.Spp)),colormap(cmocean('thermal')),caxis([-2 1])
ylabel('$f$ [Hz]','interpreter','latex')
cb = colorbar;
ylabel(cb,'[m$^2$/Hz]','interpreter','latex')
ax = gca;
set(ax,'ydir','normal','ticklabelinterpreter','latex','tickdir','out','plotboxaspectratio',[1.5 1 1])
figname = sprintf('%s/figures/%s_spectra.pdf',outRoot,L1FRoot);
exportgraphics(fig,figname)


