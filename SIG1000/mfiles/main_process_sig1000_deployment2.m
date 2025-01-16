clear all
close all
% stages of processing
% 1) define deployment number:
deploy  = 2;
% need to shift time to convert EST to UTC
time_shift = 5/24;
% 2) raw data input directory & filename convention:
rootDIR = sprintf('/Users/derekgrimes/OneDriveUNCW/DATA/BOEM/FPSD%d_SIG1000/mat_data/',deploy);
fRoot   = sprintf('S103071A011_FPS%d_',deploy);
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
atmosphTime = [datenum('15-Feb-2024 13:12:05'), datenum('15-Feb-2024 14:08:32')];
% 6) deploy/recovery times
deployTime  = [datenum('15-Feb-2024 15:00:00')]; %datenum('09-Oct-2023 16:00:00');
recoverTime = [datenum('17-Mar-2024 15:00:00')]; %datenum('30-Oct-2023 14:00:00');
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
% $$$ load_and_process_sig1000_RDI_matrix_format
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


%
Time = datetime(currents.Time,'convertFrom','datenum');
P = currents.Pressure;
U = currents.Velocity_East;
V = currents.Velocity_North;
W = currents.Velocity_Up;
qc= currents.qcFlag;
%
U(qc<0.95)=nan;
V(qc<0.95)=nan;
W(qc<0.95)=nan;
%
fig = figure,
ax1 = subplot(3,1,1);
imagesc(Time,currents.bin_mab,U,'alphadata',~isnan(U)),colormap(cmocean('balance'))
hold on, plot(Time,P,'-k','linewidth',2)
ylabel(ax1,'$(z+h)$ [m]')
caxis(ax1,[-0.5 0.5]),
cb1 = colorbar;
ylabel(cb1,'$u$ [m/s]','interpreter','latex')
set(ax1,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'ylim',[0 max(P)])
%
ax2 = subplot(3,1,2);
imagesc(Time,currents.bin_mab,V,'alphadata',~isnan(V)),colormap(cmocean('balance'))
hold on, plot(Time,P,'-k','linewidth',2)
ylabel('$(z+h)$ [m]')
caxis(ax2,[-0.25 0.25]),
cb2 = colorbar;
ylabel(cb2,'$v$ [m/s]','interpreter','latex')
set(ax2,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','xticklabel',[],'ylim',[0 max(P)])
%
ax3 = subplot(3,1,3);
imagesc(Time,currents.bin_mab,W,'alphadata',~isnan(W)),colormap(cmocean('balance'))
hold on, plot(Time,P,'-k','linewidth',2)
ylabel('$(z+h)$ [m]')
caxis(ax3,[-0.05 0.05]),
cb3 = colorbar;
ylabel(cb3,'$w$ [m/s]','interpreter','latex')
set(ax3,'tickdir','out','ticklabelinterpreter','latex','ydir','normal','ylim',[0 max(P)])
%
figname = sprintf('%s/figures/%s_currents_depth_varying.pdf',outRoot,L1FRoot);
exportgraphics(fig,figname)


Ubar = nanmean(U,1);
Vbar = nanmean(V,1);
Wbar = nanmean(W,1);

fig = figure;
plot(Time, Ubar,'-r',Time,Vbar,'-b',Time,Wbar,'-k','linewidth',2)
ylabel(' [m/s]','interpreter','latex')
h = legend('$\bar{u}$','$\bar{v}$','$\bar{w}$');
set(h,'interpreter','latex','orientation','horizontal')
set(gca,'tickdir','out','ticklabelinterpreter','latex','plotboxaspectratio',[1 1/2 1],'xlim',[Time(1) Time(end)])
figname = sprintf('%s/figures/%s_currents_depth_averaged.pdf',outRoot,L1FRoot);
exportgraphics(fig,figname)
