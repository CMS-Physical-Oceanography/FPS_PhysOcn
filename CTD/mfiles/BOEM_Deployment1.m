% "rosette", "cape-fear", "both", "seahawk", "flow-through" or "all"?
ctdType = 'seahawk';
rootDir = '/Users/derekgrimes/OneDriveUNCW/Documents-UNCW-BOEM-FryingPanShoals/General/data/BOEM_deployment1/';
rawDir  = [rootDir, filesep, 'RVSH_flowthrough',filesep];
arcDir  = [rawDir , filesep, 'L0',filesep];
if ~exist(arcDir, 'dir')
    mkdir(arcDir)
end
%
% if there are subdirectories based on dates, 
dates={''};
data = convert_seaHawk_flowThrough_to_mat(rawDir,dates,arcDir);


time = datetime(data(1).time,'convertFrom','datenum');
figure, plot(time, data(1).fDOM )
ylabel('FDOM [QSU]','interpreter','latex')
set(gca,'tickdir','out','ticklabelinterpreter','latex')
figname=sprintf('Deploy1_RVSH_flowthrough_%s_FDOM.pdf',datestr(data(1).time(1),'yyyymmdd'));
exportgraphics(gcf,figname)

tlim = time(1)+[0 0.05];
set(gca,'xlim', tlim)

figname=sprintf('Deploy1_RVSH_flowthrough_%s_FDOM_zoom.pdf',datestr(data(1).time(1),'yyyymmdd'));
exportgraphics(gcf,figname)
