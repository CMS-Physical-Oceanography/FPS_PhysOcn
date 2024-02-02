%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a 6-panel plot of:
% p  t
% u  s
% v  chl
% steps include:
% 1) path to each data-set: FPSE1, FPSS1, FPSC1
% 2) load specific fields for each mooring
%    - attempt to rescale pressure at S1?
%    - depth average, estimate %var explained
% 3) figure size/position parameters
% 4) plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
% 1) paths
rootDir = '/Users/derekgrimes/OneDriveUNCW/Documents-UNCW-BOEM-FryingPanShoals/General/data/BOEM_deployment1/';
e1dir = [rootDir,'FPSE1/L0/'];
s1dir = [rootDir,'FPSS1/L0/'];
c0dir = [rootDir,'FPSC0/L0/'];
% 2) data from each instrument
e1opt = load([e1dir,'RBR_00214016_DEP1_FPSE1_L0.mat']);
e1vel = load([e1dir,'RDI_00006104_DEP1_FPSE1_L0.mat']);
e1ctd = load([e1dir,'SBE_00003571_DEP1_FPSE1_L0.mat']);
%
s1opt = load([s1dir,'RBR_00214015_DEP1_FPSS1_L0.mat']);
s1vel = load([s1dir,'RDI_00011548_DEP1_FPSS1_L0.mat']);
s1ctd = load([s1dir,'SBE_00001308_DEP1_FPSS1_L0.mat']);
%
c0opt = load([c0dir,'0C6_23600142_DEP1_FPSC0_L0.mat']);
c0vel = load([c0dir,'SIG_00103071_DEP1_FPSC0_L0.mat']);
c0ctd = load([c0dir,'SBE_00003570_DEP1_FPSC0_L0.mat']);
%
% 2a) rescale e1vel.currents.pressure
e1_time = e1vel.currents.time';
e1_pres = e1vel.currents.pressure';
s1_pres = interp1(s1vel.currents.time,s1vel.currents.pressure,e1_time);
c0_pres = interp1(c0vel.Time         ,c0vel.Pressure         ,e1_time);
avgPres = (s1_pres+c0_pres)/2;
% least squares fit independent vars
L2var    = [0*e1_time+1, e1_pres];
% remove nans
badInds = isnan(e1_pres) | isnan(avgPres);
L2var   (badInds,:)=[];
avgPres(badInds,:)=[];
% do fit
L2fit   = L2var\avgPres;
% rescale pressure using fit
e1vel.currents.pressure = L2var(1) + L2var(2)*e1vel.currents.pressure;
%
% 2b) extract currents, apply qc flag, and depth avg
e1u = e1vel.currents.velocity_east; e1qc =~e1vel.currents.qcFlag; e1U = nansum(e1u.*e1qc,1)./sum(e1qc,1);
s1u = s1vel.currents.velocity_east; s1qc =~s1vel.currents.qcFlag; s1U = nansum(s1u.*s1qc,1)./sum(s1qc,1);
c0u = c0vel.VelEast               ; c0qc = c0vel.qcFlag         ; c0U = nansum(c0u.*c0qc,2)./sum(c0qc,2);
%
e1v = e1vel.currents.velocity_east;  e1V = nansum(e1v.*e1qc,1)./sum(e1qc,1);
s1v = s1vel.currents.velocity_east;  s1V = nansum(s1v.*s1qc,1)./sum(s1qc,1);
c0v = c0vel.VelNorth              ;  c0V = nansum(c0v.*c0qc,2)./sum(c0qc,2);
%
% 2c) estimate fraction of variance explained by barotropic vel...
e1var = nansum(nansum( (e1u.^2 + e1v.^2).*e1qc ))./nansum(nansum(e1qc));
s1var = nansum(nansum( (s1u.^2 + s1v.^2).*s1qc ))./nansum(nansum(s1qc));
c0var = nansum(nansum( (c0u.^2 + c0v.^2).*c0qc ))./nansum(nansum(c0qc));
%
e1R2  = 1-nansum(nansum( ((e1u-e1U).^2 + (e1v-e1V).^2).*e1qc ))./(nansum(nansum( e1qc ))*e1var)
s1R2  = 1-nansum(nansum( ((s1u-s1U).^2 + (s1v-s1V).^2).*s1qc ))./(nansum(nansum( s1qc ))*s1var)
c0R2  = 1-nansum(nansum( ((c0u-c0U).^2 + (c0v-c0V).^2).*c0qc ))./(nansum(nansum( c0qc ))*c0var)
%
% 3) figure params (units=cm)
xm = 2; ym = 2; pw = 5; ph = 1.75; ag = 0.33;
ppos11 = [xm            ym+2*ph+2*ag pw ph];
ppos21 = [xm            ym+1*ph+1*ag pw ph];
ppos31 = [xm            ym           pw ph];
ppos12 = [xm+pw+xm ym+2*ph+2*ag pw ph];
ppos22 = [xm+pw+xm ym+1*ph+1*ag pw ph];
ppos32 = [xm+pw+xm ym           pw ph];
ps     = [2.5*xm+2*pw  1.5*ym+3*ph+2*ag];
%
tick_locations = [datenum('10-Oct-2023 00:00:00'), ...
                  datenum('15-Oct-2023 00:00:00'),...
                  datenum('20-Oct-2023 00:00:00'),...
                  datenum('25-Oct-2023 00:00:00'),...
                  datenum('30-Oct-2023 00:00:00')];
time_limits   = tick_locations([1 end]);
%
% 4) make figure and plot
fig = figure('units','centimeters','color','w');
pos = get(fig,'position');
pos(3:4)=ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'inverthardcopy','off')
%
a11 = axes('units','centimeters','position',ppos11);
plot( s1vel.currents.time, s1vel.currents.pressure-nanmean(s1vel.currents.pressure), '-k',...
      c0vel.Time         , c0vel.Pressure-nanmean(c0vel.Pressure)                  , '-b',...      e1vel.currents.time, e1vel.currents.pressure-nanmean(e1vel.currents.pressure), '-r',...
      'linewidth',1.5)
grid on
ylabel('$\eta$ [m]','interpreter','latex')
set(a11,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
       'xlim',time_limits,'xtick',tick_locations)
%
a21 = axes('units','centimeters','position',ppos21);
plot( s1vel.currents.time, s1U, '-k',...
      c0vel.Time         , c0U, '-b',...
      e1vel.currents.time, e1U, '-r',...
      'linewidth',1.5)
grid on
ylabel('$\overline{U}$ [m]','interpreter','latex')
set(a21,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
       'xlim',time_limits,'xtick',tick_locations)
%
a31 = axes('units','centimeters','position',ppos31);
plot( s1vel.currents.time, s1V, '-k',...
      c0vel.Time         , c0V, '-b',...
      e1vel.currents.time, e1V, '-r',...
      'linewidth',1.5)
grid on
ylabel('$\overline{V}$ [m]','interpreter','latex')
set(a31,'tickdir','out','ticklabelinterpreter','latex',...
        'xlim',time_limits,'xtick',tick_locations)
datetick('x','mm/dd','keepticks','keeplimits'),
%%%
a12 = axes('units','centimeters','position',ppos12);
plot( s1vel.currents.time, s1vel.currents.temperature, '-k',...
      c0vel.Time         , c0vel.Temperature, '-b',...
      e1ctd.time         , e1ctd.temperature, '-r',...
      'linewidth',1.5)
grid on
ylabel('$T$ [C]','interpreter','latex')
set(a12,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
       'xlim',time_limits,'xtick',tick_locations)
%
a22 = axes('units','centimeters','position',ppos22);
plot( s1ctd.time, s1ctd.salinity, '-k',...
      c0ctd.time, c0ctd.salinity, '-b',...      e1ctd.time, e1ctd.salinity, '-r',...
      'linewidth',1.5)
grid on
ylabel('$S$ [PSU]','interpreter','latex')
set(a22,'tickdir','out','ticklabelinterpreter','latex','xticklabel',[],...
       'xlim',time_limits,'xtick',tick_locations)
%
a32 = axes('units','centimeters','position',ppos32);
plot( s1opt.time, s1opt.chlorophyll_a/1e3, '-k',...
      c0opt.time, c0opt.chlorophyll_a, '-b',...
      e1opt.time, e1opt.chlorophyll_a/1e3, '-r',...
      'linewidth',1.5)
grid on
text(0,0.8,'CHl-a ','units','normalized','fontsize',12)
ylabel('[$\mu$g/L]','interpreter','latex')
set(a32,'tickdir','out','ticklabelinterpreter','latex',...
        'xlim',time_limits,'xtick',tick_locations)
datetick('x','mm/dd','keepticks','keeplimits'),

figname = [rootDir, 'BOEM_deployment1_PUV_TSchl-a.png'];
exportgraphics(fig,figname)