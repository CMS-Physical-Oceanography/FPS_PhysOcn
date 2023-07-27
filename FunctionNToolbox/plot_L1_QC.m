function [L1plots]plot_L1_QC(deploypath)
% plot_L1_QC takes the .mat file from the L1processing folder and plots the data
% from the ADCP, fluorometer, and CTD. 
%  -input: deploypath (path to the BOEMTest Folder within the FPS project folder)
%  -output: L1plots
%

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
    ylabel('Pres (dB)'); title('BOEM Test 1: Pressure, Temperature, & Salinity');
    hold off; axis tight; grid on; datetick('x','keepticks','keeplimits')
subplot(312)
    plot(SBE37.time, SBE37.temperature, 'LineWidth',1.5, 'Color','red'); hold on;
    plot(ADCP.cur.mtime, ADCP.cur.temperature, 'LineWidth',1.5, 'Color', 'blue');
    plot(RBRsoloT.time, RBRsoloT.temp, 'LineWidth',1.5, 'Color','green');
    ylabel('Temp (deg C)'); hold off; legend('SBE37', 'RDI', 'RBRSoloT');
    axis tight; grid on; datetick('x','keepticks','keeplimits')
subplot(313)
    plot(SBE37.time, SBE37.salinity, 'LineWidth',1.5, 'color', 'red');
    ylabel('Sal (PSU)'); axis tight; ylim([34.5 36]); grid on;
    datetick('x','keepticks','keeplimits')

%East North Up RDI Spectrum Files
ADCPuvwfig=figure
subplot(311)
    imagesc(squeeze(ADCP.cur.east_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:45:length(ADCP.cur.east_vel)]);
    set(gca, 'ylim',[1 45]); title('East velcoity (m/s)'); colorbar
    caxis([-0.1 0.1])

subplot(312)
    imagesc(squeeze(ADCP.cur.north_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:45:length(ADCP.cur.north_vel)]);
    set(gca, 'ylim',[1 45]); title('North velcoity (m/s)');  
    colorbar; caxis([-0.1 0.1])

subplot(313)
    imagesc(squeeze(ADCP.cur.vert_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:45:length(ADCP.cur.vert_vel)]);
    set(gca, 'ylim',[1 45]); title('Vertical velcoity (m/s)');  
    colorbar; caxis([-0.05 0.05])

%RBRtri data
rbrtrifig=figure(); clf;
 subplot(311);
    plot(RBRtri.time,RBRtri.chlorophyll_a, 'LineWidth',1.5,'Color','blue');
    title('RBR Tridente'); ylabel('phyll-a (ug/l)'); grid on;
    set(gca,'XTickLabel',[]); axis tight;
subplot(312);
    plot(RBRtri.time, RBRtri.FDOM, 'LineWidth',1.5,'Color','green');
    ylabel('FDOM (ppb)'); grid on; set(gca,'XTickLabel',[]);
    axis tight;
subplot(313);
    plot(RBRtri.time, RBRtri.turbidity, 'LineWidth',1.5,'Color','red');
    ylabel('turb (FTU)'); grid on; datetick('x','keepticks','keeplimits');
    axis tight;

%RDI WH Wave Data
ADCPwavefig=figure
subplot(311)
    plot(ADCP.log9.wave_time, ADCP.log9.Hs, 'Color','blue','LineWidth',1.5);
    title('RDI WH Wave Data'); ylabel('Hs (m)'); datetick('x', 'keepticks','keeplimits');
    grid on; axis tight;
subplot(312)
    plot(ADCP.log9.wave_time,ADCP.log9.Tp, 'Color','green','LineWidth',1.5);
    ylabel('Tp (s)'); datetick('x', 'keepticks','keeplimits'); grid on;
    axis tight;
subplot(313)
    plot(ADCP.log9.wave_time,ADCP.log9.Md, 'Color','black','LineWidth',1.5);
    ylabel('Md (deg N)'); datetick('x', 'keepticks','keeplimits'); grid on;
    axis tight;

%RDI Spec Files
ADCPspecfig=figure
subplot(311)
    imagesc(squeeze(ADCP.spec.PSpec.burst)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:6:length(ADCP.spec.time)]);
    set(gca, 'ylim',[1 45]); title('Pressure Spectrum'); colorbar
    clim([-25 25]);

subplot(312)
    imagesc(squeeze(ADCP.spec.SSpec.burst)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:6:length(ADCP.spec.time)]);
    set(gca, 'ylim',[1 45]); title('Surface Spectrum');  
    colorbar; clim([-25 25]);

subplot(313)
    imagesc(squeeze(ADCP.spec.VSpec.burst)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:6:length(ADCP.spec.time)]);
    set(gca, 'ylim',[1 45]); title('Velocity Spectrum');  
    colorbar; clim([-25 25]);

%Directional Spectrum where Hs is max
HsMax_ind=find(ADCP.log9.Hs == max(ADCP.log9.Hs));
Dspec_HsMax=ADCP.spec.DSpec.burst(: ,: , HsMax_ind);

ADCPdirspecfig=figure
imagesc(Dspec_HsMax); set(gca, 'Ydir','normal'); ylim([1 64]);
ylabel('Frequency Bins'); title('Directional Spectrum');
xlabel('Directional Bins'); 

% Saving file to L1processing folder
L1Plots= [inpath, filesep, 'L1Plots.mat'];
save(L1Plots, 'pres_temp_salfig', 'ADCPuvwfig', 'rbrtrifig', 'ADCPwavefig', 'ADCPdirspecfig');


 

