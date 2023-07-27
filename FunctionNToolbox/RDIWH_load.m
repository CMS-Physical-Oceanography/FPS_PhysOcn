 function [ADCP] = RDIWH_load(deploypath)

   %======================================================================
   % This function inputs the path of a deployment folder from Wavemon
   % processing. 
   % The spectrum files stored in the waves folder are read into the
   % script. It outputs a 64xNbursts array of fourier coefficients 
   % from each burst in the deployment, and the timestamp of each burst. The
   % significant wave height, peak wave period, and mean wave direction are
   % saved to the structure and plotted for the entirety of the
   % deployment.
   %
   % The Currents file (.PD0) is ran in multiple sections within the code. 
   % The main function needed for loading currents into the structure is 
   % rdradcp.m.It specifically reads the raw binary RDI BB/Workhorse ADCP file 
   % NAME and puts all the relevant configuration and measured data into a 
   % data structure ADCP.  
   % 
   % Detailed steps of user inputs to this function are found 
   % Workhorse Sentinel 1200 kHz Post-Processing (1).docx
   %======================================================================
     
  %======================== Spectrum Profiles ============================
  
  %Find RDIWH folder within the BOEMTest folder
    rdiwhFolder=dir([deploypath, filesep, 'RDI*']);
  % find the folder with wavesMON output inside of the RDIWH folder
    wavesFolder=dir([deploypath, filesep, rdiwhFolder(1).name, filesep, 'Waves*']);
    % add folders to inpath
    inpath=[deploypath, filesep, rdiwhFolder(1).name, filesep, wavesFolder(1).name, filesep];

    for m=1:length(wavesFolder)
     %Loading the Psepc, Sspec, and Vspec files into the script to
     %concatenate all the files into separate variables with the data
     %labeled as burst. 
       % get input filenames
       source={fullfile(inpath,'PSpec*'),fullfile(inpath,'SSpec*'), fullfile(inpath,'VSpec*')};
       ADCP(m)=struct();
       bur_names = ['P_burst';'S_burst';'V_burst'];
       time = ['P_time'; 'S_time';'V_time'];
       f = ['P_freq'; 'S_freq';'V_freq'];
       for i=1:length(source)
           files=dir(source{i});
           files={files.name};
           % initialize output arrays
           dim = size(files); 
           % initial timestamp array
           times = strings(1,dim(2));
           % initial burst array 
           bursts = zeros(64,dim(2));
      %Reads the length of each spetrum file type and slices data from headers                         
        for j=1:dim(2)
           file = files{1,j};
           fil=append(inpath,'\', file);
           strdata = fileread(fil);
           % slice data from header for Vspec, Pspec, Spec files 
           data0idx = strfind(strdata,'830078)'); 
           data = str2num(strdata(data0idx+14:end));
           burst = data;
           bursts(:,j)=burst;
           num_expression = '(\d+)';
           match(j) = regexp(file, num_expression, 'match');
           times(j)=datetime(match{j},'InputFormat','yyMMddHHmm');
         
        end
          %Reading files for frequency bands and units
           fid=fopen(fil); 
           dum=fgetl(fid);
           tline=fgetl(fid);
           units = regexp(tline,'(?<=units of\s).*','match');
          %Andvance to the third line of the .txt files to read freq bands
           tline=fgetl(fid);
           f_incr = str2double(regexp(tline,'(?<=are\s).*(?=\sHz)','match'));
           f_unit=regexp(tline,'(?<=0.01562500\s).*(?=\swide)','match');
           f0= str2double(regexp(tline,'(?<=at\s).*(?=\))','match'));
           num_freq= str2double(regexp(tline,'(?<=%\s).*(?=\sFrequency)','match'));
           freq=f0:f_incr:(num_freq)*f_incr;
           
       %Savings burst and time to a variable before being overwritten by
       %the next loop.
        eval(sprintf('%s = bursts;',bur_names(i,:)))
        eval(sprintf('%s = times;',time(i,:)))
        eval(sprintf('%s = freq;',f(i,:)))
        clear bursts; clear times; clear freq;
        disp('done')
       end
       
       %DSpec Files
       source=fullfile(inpath, 'DSpec*');
        files=dir(source);
           files={files.name};
           fil=append(inpath,'\', files);
           % initialize output arrays
           dim = size(files); 
           % initial timestamp array
           times = strings(1,dim(2));
           % initial burst array 
           bursts = zeros(64,90,dim(2));
          %Reading files for spectrum units
           fid=fopen(fil{1}); 
           dum=fgetl(fid); dum=fgetl(fid); dum=fgetl(fid);
           tline=fgetl(fid);
           dir_units = regexp(tline,'(?<=are\s).*(?=\per cycle)','match');
                  
        for h=1:dim(2)
           file = files{1,h};
           fil=append(inpath,'\', file);
           strdata = fileread(fil);
           dir0idx = strfind(strdata,'begins at'); % find index 
           dir0 = strdata(dir0idx+10:dir0idx+12); % slice first bin edge
           dbinedges= str2num(dir0)+[0:89]*4;
           dbinedges(dbinedges>360)= dbinedges(dbinedges>360)-360;
           
           % slice data from header 
           data0idx = strfind(strdata,'degrees'); 
           data = str2num(strdata(data0idx+10:end));
           burst = data;
           bursts(:,:,h) = burst;
           %converting file name to datetime expression
           num_expression = '(\d+)';
           match(h) = regexp(file, num_expression, 'match');
           times(h)=datetime(match{h},'InputFormat','yyMMddHHmm');
        end
         %Savings burst and time to a variable.
           eval(sprintf('%s = bursts;', 'D_burst'));
           eval(sprintf('%s = times;', 'D_time'));
    end


  disp('Spectrum Files have been concatenated into variables')
 %========================================================================
 %%
 %====================== LOG 9 TXT file Processing =======================
 % The waves data stored within the log9 file is read in, and th data is
 % stored in the ADCP structure under the field log9
  for mm=1:length(wavesFolder)  
    logfile_source=fullfile(inpath,'*LOG9.TXT');
    logfile_dir=dir(logfile_source);
    file=append(inpath,'\', logfile_dir.name);
    fileID = fopen(file);
    [C] = textscan(fileID, '%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f , %f, %f, %f, %f, %f, %f , %f, %f, %f %*[^\n]');  
    fclose(fileID);
   %Creating time vector
    year = C{1,2};
    month = C{1,3};
    day = C{1,4};
    hour = C{1,5};
    minu = C{1,6};
    sec = C{1,7};
    wave_time= datetime([year+2000 month day hour minu sec]); %adding 2000 to make YY into YYYY format

    H_s = C{1,9}; %Significant Wave Height (meters)
    T_p = C{1,10}; %Peak Wave Period (seconds) - period associated with the largest peak in the power spectrum
    D_p = C{1,11}; %Peak Wave Direction (degrees) - peak direction at the peak period.
    Tp_Sea = C{1,12}; %Peak Sea Wave Period (seconds) - period associated with the largest peak in the sea region of the power spectrum.
    Dp_Sea = C{1,13}; % Peak Sea Wave Direction (degrees) - peak sea direction at the peak period in the sea region.
    Hs_Sea = C{1,14}; % Significant Wave Height in the sea region of the power spectrum
    Tp_Swell = C{1,15}; % Peak Swell Wave Period (seconds) - period associated with the largest peak in the swell region of the power spectrum.
    Dp_Swell = C{1,16}; % Peak Swell Wave Direction (degrees) - peak swell direction at the peak period in the swell region.
    Hs_Swell = C{1,17}; % Significant Wave Height in the swell region of the power spectrum
    logDepth = C{1,18}; % Water level (from pressure sensor) (millimeters)
    Hmax = C{1,19}; % Maximum wave height (meters) as determined by Zero-Crossing analysis of the surface track time series. 
    Tmax = C{1,20}; % Maximum Peak Wave Period (seconds) as determined by Zero-Crossing analysis of the surface track time series
    Honethird = C{1,21}; % Significant wave height of the largest 1/3 of the waves in the field as determined by Zero-Crossing analysis of the surface track time series.
    Tonethird = C{1,22}; % The period associated with the peak wave height of the largest 1/3 of the waves in the field as determined by Zero-Crossing analysis of the surface track time series.
    Hmean = C{1,23}; % The mean significant wave height of the waves in the field as determined by ZeroCrossing analysis of the surface track time series
    Tmean = C{1,24}; % The period associated with the mean significant wave height of the waves in the field as determined by Zero-Crossing analysis of the surface track time series
    Honetenth = C{1,25}; % Significant wave height of the largest 1/10 of the waves in the field as determined by Zero-Crossing analysis of the surface track time series.
    Tonetenth = C{1,26}; % The period associated with the peak wave height of the largest 1/10 of the waves in the field as determined by Zero-Crossing analysis of the surface track time series.
    Dmean = C{1,27}; % Mean Peak Wave Direction (degrees)

    
    %Commented out. Needed if the user wants to modify the raw data within
    %the function.
    %Plotting Raw data to chop out of water time
%     figure(); clf
%     plot(datenum(wave_time), H_s, 'LineWidth', 1.5,'Color','blue'); grid on;
%     ylabel('Hs (m)'); title('Select a Start & End Time with Cursor');
    
   
%     [Hs_xx,Hs_yy]=ginput(2);
%     
%     Hs_predata=find(datenum(wave_time) < Hs_xx(1));
%     Hs_postdata=find(datenum(wave_time) > Hs_xx(2));
%     
%   %Post Deployment data delete
%    H_s(Hs_postdata)=[];
%    T_p(Hs_postdata)=[];
%    M_d(Hs_postdata)=[];
%    wave_time(Hs_postdata)=[];
%   %Pre-Deployment data delete
%    H_s(Hs_predata)=[];
%    T_p(Hs_predata)=[];
%    M_d(Hs_predata)=[]; 
%    wave_time(Hs_predata)=[];
   
    %Plotting Wave Data Time Series
    figure(); clf
    subplot(311)
        plot(wave_time, H_s,'LineWidth', 1.5,'Color','blue'); grid on
        set(gca, 'XTickLabel',[]); ylabel('Hs (m)'); title('Raw Wave Data'); 
        axis tight; 
    subplot(312)
        plot(wave_time, T_p,'LineWidth', 1.5,'Color','green'); grid on
        set(gca, 'XTickLabel',[]); ylabel('Tp (s)'); axis tight;
    subplot(313)
        plot(wave_time, D_p, '.','LineWidth', 1.5,'Color','magenta'); 
        grid on; ylabel('D_p (deg)'); axis tight; datetick('x','dd/mm','keeplimits','keepticks')
  end
        
%% Saving Wave and Spectrum data to ADCP structure

  %Data Structure stroing the Spectrum samples 
     ADCP.spec.PSpec.burst=P_burst;
     ADCP.spec.PSpec.units=units;
     ADCP.spec.PSpec.freq=P_freq';
     ADCP.spec.SSpec.burst=S_burst;
     ADCP.spec.SSpec.units=units;
     ADCP.spec.SSpec.freq=S_freq';
     ADCP.spec.VSpec.burst=V_burst;
     ADCP.spec.VSpec.units=units;
     ADCP.spec.VSpec.freq=V_freq';
     ADCP.spec.DSpec.burst=D_burst;
     ADCP.spec.DSpec.units=dir_units;
     ADCP.spec.DSpec.dir_bins=dbinedges;
     ADCP.spec.DSpec.freq=P_freq';
     ADCP.spec.time=convertTo(datetime(P_time, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'), 'datenum')';
 % Data Structure for Log 9 File
     ADCP.log9.Hs=H_s;
     ADCP.log9.Tp=T_p;
     ADCP.log9.Dp=D_p;
     ADCP.log9.Tp_Sea = Tp_Sea;
     ADCP.log9.Dp_Sea = Dp_Sea;
     ADCP.log9.Hs_Sea = Hs_Sea;
     ADCP.log9.Tp_Swell = Tp_Swell;
     ADCP.log9.Dp_Swell = Dp_Swell;
     ADCP.log9.Hs_Swell = Hs_Swell;
     ADCP.log9.Hmax = Hmax;
     ADCP.log9.Tmax = Tmax;
     ADCP.log9.Honethird = Honethird;
     ADCP.log9.Tonethird = Tonethird;
     ADCP.log9.Honetenth = Honetenth;
     ADCP.log9.Tonetenth = Tonetenth;
     ADCP.log9.Hmean = Hmean;
     ADCP.log9.Tmean = Tmean;
     ADCP.log9.D_pmean=Dmean
     ADCP.log9.units={'Hs= m'; 'Tp = s'; 'Md = deg'};
     ADCP.log9.wave_time=convertTo(wave_time, 'datenum');

%%       
%============================= Currents =================================
% Using the rdradcp to go from an ADCP binary to a mat file from a .PD0 file
% If ADCP was recordeing waves you must process the raw files with waves-
% mon first, then use the resulting .PD0 file here. Code Modified from Cody Brenton

 
%Section 1: 
  wave_source=fullfile(inpath,'*.PD0');
  pd0file=dir(wave_source);
  cd(pd0file.folder);
  pd0filedir=fullfile(pd0file.folder, pd0file.name)

  ct = 0;
  for h = 1:length(pd0file) 
    try
        [ADCP.cur(h)] = rdradcp(pd0file(h).name,1,-1)
    catch
        try
            [ADCP.cur(h)] = rdradcp(pd0file(h).name,1,[1 16000])
        catch
            try
               [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 15000])
            catch
                try
                    [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 14000])
                catch
                    try
                        [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 13000])
                    catch
                        try
                            [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 12000])
                        catch
                            try
                                [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 11000])
                            catch
                                try
                                    [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 10000])
                                catch
                                    try
                                        [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 9000])
                                    catch
                                        try
                                            [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 8500])
                                        catch
                                            try
                                                [ADCP.cur(h)]=rdradcp(pd0file(h).name,1,[1 8000])
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    ct = ct+1
  end
  disp('Section 1 complete')
%% Section 2 
% The bad values will have an mtime = 0 from the previous loop assinged
% from the rdradcp function.

ct=0;
for i = 1:length(ADCP.cur)

    index = find(ADCP.cur(i).mtime == 0);
    ADCP.cur(i).mtime(index) = [];
    ADCP.cur(i).number(index)=[];
    ADCP.cur(i).pitch(index)=[];
    ADCP.cur(i).roll(index)=[];
    ADCP.cur(i).heading(index)=[];
    ADCP.cur(i).pitch_std(index)=[];
    ADCP.cur(i).roll_std(index)=[];
    ADCP.cur(i).heading_std(index)=[];
    ADCP.cur(i).depth(index)=[];
    ADCP.cur(i).temperature(index)=[];
    ADCP.cur(i).salinity(index)=[];
    ADCP.cur(i).pressure(index)=[];
    ADCP.cur(i).pressure_std(index)=[];
    ADCP.cur(i).east_vel(:,index)=[];
    ADCP.cur(i).north_vel(:,index)=[];
    ADCP.cur(i).vert_vel(:,index)=[];
    ADCP.cur(i).error_vel(:,index)=[];
    ADCP.cur(i).corr(:,:,index)=[];
    ADCP.cur(i).status(:,:,index)=[];
    ADCP.cur(i).intens(:,:,index)=[];
    ADCP.cur(i).bt_range(:,index)=[];
    ADCP.cur(i).bt_vel(:,index)=[];
    ADCP.cur(i).bt_corr(:,index)=[];
    ADCP.cur(i).bt_ampl(:,index)=[];
    ADCP.cur(i).bt_perc_good(:,index)=[];
    
    ct=ct+1;
end
disp('Section 2 complete')
%%  Section 3 Test Images 
% creates images to visualize the data for the first and last 100 ensembles.
figure
for jj=1:length(ADCP.cur)
    subplot(221)
    imagesc(squeeze(ADCP.cur(jj).intens(:,1,:)));
    title('Beam 1 Echo Amp'); ylabel('bin number');
    set(gca, 'Ydir', 'normal'); xticks([1:75:length(ADCP.cur.intens)]); grid on;
    
    ax2=subplot(222)
    imagesc(squeeze(ADCP.cur(jj).intens(:,2,:)));
    title('Beam 2 Echo Amp');ylabel('bin number');
    set(gca, 'Ydir', 'normal'); xticks([1:75:length(ADCP.cur.intens)]); grid on;
    colorbar(ax2); caxis([40 160]);

    subplot(223)
    imagesc(squeeze(ADCP.cur(jj).intens(:,3, :)));
    title('Beam 3 Echo Amp');ylabel('bin number');
    set(gca, 'Ydir', 'normal'); xticks([1:75:length(ADCP.cur.intens)]); grid on;

    ax4=subplot(224)
    imagesc(squeeze(ADCP.cur(jj).intens(:,4, :)));
    title('Beam 4 Echo Amp');ylabel('bin number'); set(gca, 'ylim',[1 45]);
    set(gca, 'Ydir', 'normal'); xticks([1:75:length(ADCP.cur.intens)]); grid on;
    colorbar(ax4); caxis([40 160]);
end

figure 
for kk=1:length(ADCP.cur)
    subplot(221)
    imagesc(squeeze(ADCP.cur(kk).east_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:75:length(ADCP.cur.east_vel)]);
    set(gca, 'ylim',[1 45]); title('East velcoity (m/s)'); colorbar
    caxis([-0.25 0.25])

    subplot(222)
    imagesc(squeeze(ADCP.cur(kk).north_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:75:length(ADCP.cur.north_vel)]);
    set(gca, 'ylim',[1 45]); title('North velcoity (m/s)');  
    colorbar; caxis([-0.25 0.25])

    subplot(223)
    imagesc(squeeze(ADCP.cur(kk).vert_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:75:length(ADCP.cur.vert_vel)]);
    set(gca, 'ylim',[1 45]); title('Vertical velcoity (m/s)');  
    colorbar; caxis([-0.15 0.15])

    subplot(224)
    imagesc(squeeze(ADCP.cur(kk).error_vel)); set(gca, 'YDir','normal');
    ylabel('bin number'); xticks([1:75:length(ADCP.cur.error_vel)]);
    set(gca, 'ylim',[1 45]); title('Error velcoity (cm/s)');  
    colorbar; caxis([-0.1 0.1])
end

figure
for hh=1:length(ADCP.cur)
    subplot(311)
    plot(ADCP.cur(hh).mtime,ADCP.cur(hh).heading, 'Color','blue','LineWidth',1.5);
    title('Raw Directional Pivoting during Deployment'); ylabel('Heading')

    subplot(312)
    plot(ADCP.cur(hh).mtime,ADCP.cur(hh).pitch, 'Color','green','LineWidth',1.5);
    ylabel('Pitch');

    subplot(313)
    plot(ADCP.cur(hh).mtime,ADCP.cur(hh).roll, 'Color','black','LineWidth',1.5);
    ylabel('Roll');
end
disp('Section 3 complete')

%% Section 4 Calculating depth from pressure data
%the provided depth values are discrete so here we calculate our own depth
%field and call it "new_depth"
for ii = 1:length(ADCP.cur)
    if isempty(ADCP.cur(ii).pressure) == 0
        for jj = 1:length(ADCP.cur(ii).mtime)
            ADCP.cur(ii).mean_depth(jj) = nanmean((ADCP.cur(ii).pressure(jj))*10/(1026*9.81)); %recorded pressure is in decibars so we use density of water to be 1026 and g to be 9.81
        end
    else
        ADCP.cur(ii).mean_depth = []; %if there was no valid pressure data the field is empty
    end
end
disp('Section 4 complete')

%% Section 5 Remove fields that are not needed for L2

fields = {'pitch_std','roll_std','heading_std','status','bt_range','bt_vel',...
    'bt_corr','bt_ampl','bt_perc_good','perc_good', 'number', 'pressure_std','salinity', 'coords', 'rerunparameter'};

ADCP.cur = rmfield(ADCP.cur,fields);
disp('Section 5 complete')

%% Save
%Saving the file to outpath directory with file name and start date of
%deployment
% Creating L0 Folder within the BOEMTest and Crate Folder
L0Folder=mkdir(deploypath, 'L0processing')

%Creating path to the L0 folder within BoemTest Folder
outpath=[inpath,'..', filesep, '..'];
cd(outpath); outpath=pwd;
L0folder=dir([outpath, filesep,'L0*']);
outpath=[outpath, filesep, L0folder(1).name];
%Adding Date of Deployment and RDIWH SN to filename
 date=datestr(ADCP.spec.time(1),'yyyy_mm_dd');
 rdiSN=num2str(ADCP.cur.config.remus_serialnum);
RDIfilename=[outpath, filesep,'RDIWH_',rdiSN,'_', date,'L0.mat'];
save(RDIfilename, 'ADCP') 
     
disp('You have successfully loaded all data from WaveMon file into the ADCP structure')
 end
            
        
        
        
               
               