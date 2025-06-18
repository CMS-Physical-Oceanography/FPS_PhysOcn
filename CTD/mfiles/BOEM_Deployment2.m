% there were errors with the GPS stream for this experiment, requiring custom processing scripts
clear all
close all
ctdType = 'seahawk';
rootDir = '/Users/derekgrimes/OneDriveUNCW/Documents-UNCW-BOEM-FryingPanShoals/General/data/BOEM_deployment2/';
rawDir  = [rootDir, filesep, 'RVSH_flowthrough',filesep];
arcDir  = [rawDir , filesep, 'L0',filesep];
if ~exist(arcDir, 'dir')
    mkdir(arcDir)
end
%
% if there are subdirectories based on dates, 
dates={''};
data = convert_seaHawk_flowThrough_to_mat_deploy2(rawDir,dates,arcDir);


% $$$ time = datetime(data(1).time,'convertFrom','datenum');
% $$$ figure, plot(time, data(1).fDOM )
% $$$ ylabel('FDOM [QSU]','interpreter','latex')
% $$$ set(gca,'tickdir','out','ticklabelinterpreter','latex')
% $$$ figname=sprintf('Deploy1_RVSH_flowthrough_%s_FDOM.pdf',datestr(data(1).time(1),'yyyymmdd'));
% $$$ exportgraphics(gcf,figname)
% $$$ 
% $$$ tlim = time(1)+[0 0.05];
% $$$ set(gca,'xlim', tlim)
% $$$ 
% $$$ figname=sprintf('Deploy1_RVSH_flowthrough_%s_FDOM_zoom.pdf',datestr(data(1).time(1),'yyyymmdd'));
% $$$ exportgraphics(gcf,figname)


figDir = [rawDir,filesep,'figures'];
if ~exist(figDir,'dir')
    eval(['!mkdir -p ',figDir])
end

for kk = 1:length(data)
    figure,
    a1 = subplot(8,1,1);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).temperature,'.k','markersize',8)
    ylabel(a1,'[$^\circ$C]','interpreter','latex')
    set(a1,'xticklabel',[])

    a2 = subplot(8,1,2);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).salinity,'.k','markersize',8)
    ylabel(a2,'[psu]','interpreter','latex')
    set(a2,'xticklabel',[])

    a3 = subplot(8,1,3);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).turbidity,'.k','markersize',8)
    ylabel(a3,'[ftu]','interpreter','latex')
    set(a3,'xticklabel',[])

    a4 = subplot(8,1,4);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).ChlA,'.k','markersize',8)
    ylabel(a4,'[ChlA]','interpreter','latex')
    set(a4,'xticklabel',[])
    
    a5 = subplot(8,1,5);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).fDOM,'.k','markersize',8)
    ylabel(a5,'[fDOM]','interpreter','latex')
    set(a5,'xticklabel',[])

    a6 = subplot(8,1,6);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).O2mg_L,'.k','markersize',8)
    ylabel(a6,'O$_2$ [$\mu$g/L]','interpreter','latex')
    set(a6,'xticklabel',[])

    a7 = subplot(8,1,7);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).pH,'.k','markersize',8)
    ylabel(a7,'[pH]','interpreter','latex')
    set(a7,'xticklabel',[])

    a8 = subplot(8,1,8);
    plot( datetime(data(kk).time,'convertfrom','datenum'), data(kk).pressure,'.k','markersize',8)
    ylabel(a8,'[decibar]','interpreter','latex')
    orient(gcf,'tall')

    figname = [figDir,filesep,'seahawk_flowthrough_all_vs_time.pdf'];
    exportgraphics(gcf,figname)
    
    if isempty(data(kk).latitude)
        continue
    end
    
    figure,
    geoscatter(data(kk).latitude,data(kk).longitude,10,data(kk).temperature),colormap(cmocean('thermal'))
    geobasemap('satellite'), colorbar

    figname = [figDir,filesep,'seahawk_flowthrough_temp_vs_lat_lon.pdf'];
    exportgraphics(gcf,figname)


    figure,
    geoscatter(data(kk).latitude,data(kk).longitude,10,data(kk).salinity),colormap(cmocean('thermal'))
    geobasemap('satellite'), colorbar

    figname = [figDir,filesep,'seahawk_flowthrough_salt_vs_lat_lon.pdf'];
    exportgraphics(gcf,figname)

    figure,
    geoscatter(data(kk).latitude,data(kk).longitude,10,data(kk).fDOM),colormap(cmocean('thermal'))
    geobasemap('satellite'), colorbar

    figname = [figDir,filesep,'seahawk_flowthrough_fDOM_vs_lat_lon.pdf'];
    exportgraphics(gcf,figname)

end