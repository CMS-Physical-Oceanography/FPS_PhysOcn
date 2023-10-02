function [L1plots]=plot_L1_QC(deploypath)
% plot_L1_QC takes the .mat file from the L1processing folder and plots the data
% from the ADCP, fluorometer, and CTD. 
%  -input: deploypath (path to the BOEMTest Folder within the FPS project folder)
%  -output: L1plots

%Loading in L1proccesed data structures from the L1 folder
L1folder=dir([deploypath, filesep, 'L1pr*']);
inpath=[deploypath, filesep, L1folder(1).name];
L1file=dir(fullfile(inpath, 'BOEML1*'));
load([inpath, filesep,L1file.name]);

%Pressure, Temperature, and Salinity (RDI, SBE37, RBRsoloT)
pres_temp_salfig=figure
subplot(311)
    plot(SBE37.time, SBE37.pressure, 'LineWidth',1.5,'Color','red'); hold on;
    plot(ADCP.cur.mtime, ADCP.cur.pressure./1010, 'LineWidth',1.5,'Color','blue');
    ylabel('Pres (dB)'); title(sprintf('BOEM, %s: Pressure, Temperature, & Salinity',char(inpath(end-18:end-13))));
    hold off; axis tight; grid on; datetick('x','keepticks','keeplimits')
subplot(312)
    plot(SBE37.time, SBE37.temperature, 'LineWidth',1.5, 'Color','red'); hold on;
    plot(ADCP.cur.mtime, ADCP.cur.temperature, 'LineWidth',1.5, 'Color', 'blue');
    plot(RBRsoloT.time, RBRsoloT.temp, 'LineWidth',1.5, 'Color','green');
    ylabel('Temp (deg C)'); hold off;
    axis tight; grid on; datetick('x','keepticks','keeplimits')
    legend(sprintf('SBE37: %s',char(SBE37.SN)), sprintf('RDI: %d',ADCP.cur.config.remus_serialnum), sprintf('RBRSoloT: %s', RBRsoloT.SN) ...
        ,'orientation','horizontal','location','northoutside');
subplot(313)
    plot(SBE37.time, SBE37.salinity, 'LineWidth',1.5, 'color', 'red');
    ylabel('Sal (PSU)'); axis tight; ylim([34.5 36]); grid on;
    datetick('x','keepticks','keeplimits')

%East North Up & Error Vel RDI Spectrum Files
ADCPuvwfig=figure
    subplot(221)
    imagesc(ADCP.cur.mtime, ADCP.cur.config.ranges, squeeze(ADCP.cur.east_vel)); 
    set(gca, 'YDir','normal'); ylabel('Depth (m)'); title(sprintf('BOEM, %s: East Velocity (m/s)',char(inpath(end-18:end-13)))); colorbar
    caxis([-0.25 0.25]); datetick('x', 'mm/dd','keeplimits','keepticks'); ylim([0 ADCP.cur.depth_ea+2]);

    subplot(222)
    imagesc(ADCP.cur.mtime, ADCP.cur.config.ranges, squeeze(ADCP.cur.north_vel)); 
    set(gca, 'YDir','normal'); title('North velcoity (m/s)');  ylim([0 ADCP.cur.depth_ea+2]);
    colorbar; caxis([-0.25 0.25]); datetick('x', 'mm/dd','keeplimits','keepticks');

    subplot(223)
    imagesc(ADCP.cur.mtime, ADCP.cur.config.ranges, squeeze(ADCP.cur.vert_vel)); 
    set(gca, 'YDir','normal'); ylabel('Depth (m)'); title('Vertical velcoity (m/s)');  
    colorbar; caxis([-0.005 0.005]); datetick('x', 'mm/dd','keeplimits','keepticks');
    ylim([0 ADCP.cur.depth_ea+2]);

    subplot(224)
    imagesc(ADCP.cur.mtime, ADCP.cur.config.ranges, squeeze(ADCP.cur.error_vel)); 
    set(gca, 'YDir','normal'); title('Error velcoity (cm/s)');  
    colorbar; caxis([-0.01 0.01]); datetick('x', 'mm/dd','keeplimits','keepticks');
    ylim([0 ADCP.cur.depth_ea+3]);

%RBRtri data
rbrtrifig=figure(); clf;
 subplot(311);
    plot(RBRtri.time,RBRtri.chlorophyll_a, 'LineWidth',1.5,'Color','blue');
    title(sprintf('BOEM, %s: RBR Tridente',char(inpath(end-18:end-13)))); ylabel('phyll-a (ug/l)'); grid on;
    set(gca,'XTickLabel',[]); axis tight;
subplot(312);
    plot(RBRtri.time, RBRtri.FDOM, 'LineWidth',1.5,'Color','green');
    ylabel('FDOM (ppb)'); grid on; set(gca,'XTickLabel',[]);
    axis tight;
subplot(313);
    plot(RBRtri.time, RBRtri.turbidity, 'LineWidth',1.5,'Color','red');
    ylabel('turb (FTU)'); grid on; datetick('x','mm/dd','keepticks','keeplimits');
    axis tight;

%RDI WH Wave Data
ADCPwavefig=figure
subplot(311)
    plot(ADCP.log9.wave_time, ADCP.log9.Hs, 'Color','blue','LineWidth',1.5);
    title(sprintf('BOEM, %s: RDI WH Wave Data',char(inpath(end-18:end-13)))); ylabel('Hs (m)'); datetick('x', 'keepticks','keeplimits');
    grid on; axis tight;
subplot(312)
    plot(ADCP.log9.wave_time,ADCP.log9.Tp, 'Color','green','LineWidth',1.5);
    ylabel('Tp (s)'); datetick('x', 'keepticks','keeplimits'); grid on;
    axis tight;
subplot(313)
    plot(ADCP.log9.wave_time,ADCP.log9.Dp, 'Color','black','LineWidth',1.5);
    ylabel('Md (deg N)'); datetick('x', 'keepticks','keeplimits'); grid on;
    axis tight;

%RDI Spec Files
ADCPspecfig=figure
subplot(311)
    imagesc(ADCP.spec.time, ADCP.spec.PSpec.freq,squeeze(ADCP.spec.PSpec.burst)); set(gca, 'YDir','normal');
    ylabel(sprintf('BOEM, %s', char(ADCP.spec.PSpec.units))); datetick('x','mm/dd', 'keeplimits','keepticks');
    title(sprintf('%s: Pressure Spectrum',char(inpath(end-18:end-13)))); colorbar; clim([-5 5]);

subplot(312)
    imagesc(ADCP.spec.time, ADCP.spec.SSpec.freq,squeeze(ADCP.spec.SSpec.burst)); 
     ylabel(sprintf('%s',char(ADCP.spec.SSpec.units))); datetick('x','mm/dd', 'keeplimits','keepticks');
     title('Surface Spectrum');  set(gca, 'YDir','normal');
    colorbar; clim([-5 5]);

subplot(313)
    imagesc(ADCP.spec.time, ADCP.spec.VSpec.freq, squeeze(ADCP.spec.VSpec.burst));
    ylabel(sprintf('%s',char(ADCP.spec.VSpec.units))); datetick('x','mm/dd', 'keeplimits','keepticks');
    title('Velocity Spectrum'); set(gca, 'YDir','normal');  
    colorbar; clim([-5 5]);

%Directional Spectrum where Hs is max
HsMax_ind=find(ADCP.log9.Hs == max(ADCP.log9.Hs));
Dspec_HsMax=ADCP.spec.DSpec.burst(: ,: , HsMax_ind);

ADCPdirspecfig=figure
imagesc(ADCP.spec.DSpec.freq,ADCP.spec.DSpec.dir_bins,Dspec_HsMax); set(gca, 'Ydir','normal'); 
 title(sprintf('BOEM, %s: Directional Spectrum',char(inpath(end-18:end-13)))); colorbar; xlabel(sprintf('Frequency: %s', char(ADCP.spec.DSpec.units)));
 ylabel('Direction Bins (deg)');

% Saving file to L1processing folder
L1Plots= [inpath, filesep, 'L1Plots.mat'];
save(L1Plots, 'pres_temp_salfig', 'ADCPuvwfig', 'rbrtrifig', 'ADCPwavefig', 'ADCPdirspecfig');


 

