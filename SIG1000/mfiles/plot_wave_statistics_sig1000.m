clear all
close all
%
dtAvg = 3600;% seconds
dtEns = 160;% seconds
%
root  = '/home/derek/projects/FRF2022/mat_data/';
%
% output file format:
sigNum = 4;% 04=McCorm; 05=Grimes
outRoot = root;
outFRoot= sprintf('FRF2022_sig1000_%02d_wave_stats.mat',sigNum);
%
fprintf(['loading: %s \n'],[outRoot,outFRoot])
sig4 = load([outRoot,outFRoot]);
%
sigNum = 5;% 04=McCorm; 05=Grimes
outRoot = root;
outFRoot= sprintf('FRF2022_sig1000_%02d_wave_stats.mat',sigNum);
%
fprintf(['loading: %s \n'],[outRoot,outFRoot])
sig5 = load([outRoot,outFRoot]);
%
xlims = [min(sig4.t) datenum('23-Sep-2022 09:00:00')];
xticks = [datenum('21-Sep-2022'):1:datenum('20-Oct-2022')];
%
cm = cmocean('balance');
%
xm = 2;
ym = 2;
pw = 6;
ph = 3;
ag = 0.5;
cw = 0.25;
%
ppos1 = [xm ym pw ph];
cbpos1= [ppos1(1:2)+[ppos1(3) 0] cw 0.7*ph];
ppos2 = [xm ym+ph+ag pw ph];
cbpos2= [ppos2(1:2)+[ppos2(3) 0] cw 0.7*ph];
ppos3 = [xm ym+2*ph+2*ag pw ph];
cbpos3= [ppos3(1:2)+[ppos3(3) 0] cw 0.7*ph];
ppos4 = [xm ym+3*ph+3*ag pw ph];
cbpos4= [ppos4(1:2)+[ppos4(3) 0] cw 0.7*ph];
%
ps = [xm+pw+3*ag+cw ym+4*ag+4*ph];
%
fig = figure('color','w','units','centimeters');
pos = get(fig,'position');
pos(3:4) = ps;
set(fig,'position',pos,'papersize',ps,'paperposition',[0 0 ps],'inverthardcopy','off')
%
a4 = axes('units','centimeters','position',ppos4);
p1 = plot(sig4.t,sig4.Hs,'-k',sig5.t,sig5.Hs,':k');
hh = legend(p1,'McCormack','Grimes');
set(hh,'location','northwest','fontsize',9,'interpreter','latex')

set(a4,'xlim',xlims,'xtick',xticks)
datetick('x','mm/dd','keeplimits')
ylabel(a4,'$H_s$~[m]','interpreter','latex')
set(a4,'xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',12,'ylim',[0 1.5])
%
a40 = axes('units','centimeters','position',ppos4);
plot(sig4.t,sig4.Tm,'-b',sig5.t,sig5.Tm,':b')
ylabel(a40,'$T_m$~[s]','interpreter','latex')
set(a40,'color','none','yaxislocation','right','box','off','ycolor','b','ylim',[0 15],'xtick',[],'ticklabelinterpreter','latex')
%
%
a3 = axes('units','centimeters','position',ppos3);
plot(sig4.t,sig4.mdir(:,1),'-k',sig5.t,sig5.mdir(:,2),':k')
set(a3,'xlim',xlims,'xtick',xticks)
datetick('x','mm/dd','keeplimits')
ylabel(a3,'$\theta_m$~[$^\circ$]','interpreter','latex')
set(a3,'xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',12,'ylim',[-45 15],'ydir','normal')
%
%
a2 = axes('units','centimeters','position',ppos2);
p2 = plot(sig4.t,sig4.mspread(:,1),'-r',sig4.t,sig4.mspread(:,2),'-b',sig5.t,sig5.mspread(:,1),':r',sig5.t,sig5.mspread(:,2),':b');
set(a2,'xlim',xlims,'xtick',xticks)
datetick('x','mm/dd','keeplimits')
set(a2,'xticklabel',[],'tickdir','out','ticklabelinterpreter','latex','fontsize',12,'ylim',[0 45],'ydir','normal')
ylabel(a2,'$\sigma_\theta$ [$^\circ$]','interpreter','latex')
hh = legend([p2(1), p2(2)],'(a1,b1)','(a2,b2)');
set(hh,'location','southwest','fontsize',9,'interpreter','latex')
%
%
a1 = axes('units','centimeters','position',ppos1);
plot(sig4.t,sig4.mSxy,'-k',sig5.t,sig5.mSxy,':k')
set(a1,'xlim',xlims,'xtick',xticks)
datetick('x','mm/dd','keeplimits')
ylabel(a1,'$S_{xy}$~[m$^2\,$s$^-2$]','interpreter','latex')
xlabel('mm/dd/2013')

set(a1,'tickdir','out','ticklabelinterpreter','latex','fontsize',12,'ylim',[-0.25 0.25],'ydir','normal')

figname = sprintf('/home/derek/projects/FRF2022/figures/FRF2022_sig1000_wave_stats_trimmed.pdf',sigNum);
exportgraphics(fig,figname)

