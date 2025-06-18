 function [ADCP] = RDIWH_load_waves(inpath,save_to_path,depnum,stationID,stime,etime)
%[ADCP] = RDIWH_load_waves(inpath,save_to_path,depnum,stationID,stime,etime)
%
   %======================================================================
   % This function inputs the path of a deployment folder from Wavemon
   % processing. 
   %
   % Changed to the WAVES path that comes out of WAVESMON. Should have all
   % the .txt files for diretional spectra, SVP spectra, AND the LOG9 file for bulk statistics etc.
   %
   %
   % The spectrum files stored in the waves folder are read into the
   % script. It outputs a 64xNbursts array of fourier coefficients 
   % from each burst in the deployment, and the timestamp of each burst. The
   % significant wave height, peak wave period, and mean wave direction are
   % saved to the structure and plotted for the entirety of the
   % deployment.
   %
   % 
   % Detailed steps of user inputs to this function are found 
   % Workhorse Sentinel 1200 kHz Post-Processing (1).docx
   %======================================================================
  
  % Current directory:
%  ADir = pwd;
   
  %======================== Spectrum Profiles ============================
  
  %Find RDIWH folder within the BOEMTest folder
%    rdiwhFolder=dir([deploypath, filesep, 'RDI*']);
    rdiwhFolder=dir(inpath);

%For loop created to run the function even if the folder/file for the
%SeaBird is not included on the lander
if isempty(rdiwhFolder) == 1
disp('Can''t find the right folder with RDI WH WAVES data')
else

     %Loading the Psepc, Sspec, and Vspec files into the script to
     %concatenate all the files into separate variables with the data
     %labeled as burst. 
       % get input filenames
       specfolder = dir([inpath, '*Spec*']);
       if isempty(specfolder) == 1
       %ADCP.spec = []; % Don't do anything 
       else
 %      specfile_path = [inpath, specfolder(1).name, filesep];
 %      inpath = specfile_path;
       source={fullfile(inpath,'PSpec*'),fullfile(inpath,'SSpec*'), fullfile(inpath,'VSpec*')};
       ADCP = struct();
       bur_names = ['P_burst';'S_burst';'V_burst'];
       time = ['P_time'; 'S_time';'V_time'];
       f = ['P_freq'; 'S_freq';'V_freq'];
       for i=1:length(source)
           files=dir(source{i});
               if isempty(files) == 1
               else
           files={files.name};
           % initialize output arrays
           dim = size(files); 
           % initial timestamp array
           % times = strings(1,dim(2));
           % initial burst array 
           bursts = zeros(64,dim(2));
      %Reads the length of each spetrum file type and slices data from headers                         
            for j=1:dim(2)
               file = files{1,j};
               %file1=append(inpath,'\', file);
               file1 = [inpath file];
               strdata = fileread(file1);
               % slice data from header for Vspec, Pspec, Spec files 
               data0idx = strfind(strdata,'830078)'); 
               data = str2num(strdata(data0idx+14:end));
               burst = data;
               bursts(:,j)=burst;
               num_expression = '(\d+)';
               match(j) = regexp(file, num_expression, 'match');
               times(j)=datenum(match{j},'yymmddHHMM');

            end
          %Reading files for frequency bands and units
           fid=fopen(file1); 
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
               end
       end
disp('Done Spec')    
       %DSpec Files
       source=fullfile(inpath, 'DSpec*');
        files=dir(source);
           files={files.name};
           %fil=append(inpath,'\', files);
           fil = append(inpath, files);
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
           %fil=append(inpath,'\', file);
           fil = append(inpath, file);
           strdata = fileread(fil);
           dir0idx = strfind(strdata,'begins at'); % find index 
           dir0 = strdata(dir0idx(1)+10:dir0idx(1)+11); % slice first bin edge
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
           TIME(h) = datenum(times(h));
        end
        disp('Done DSpec')    

         %Savings burst and time to a variable.
           eval(sprintf('%s = bursts;', 'D_burst'));
           eval(sprintf('%s = times;', 'D_time'));
    %Data Structure storing the Spectrum samples 
     if exist('P_burst','var') == 1
     ADCP.spec.PSpec.burst=P_burst;
     ADCP.spec.PSpec.units=units;
     ADCP.spec.PSpec.freq=P_freq';
     else
     end
     if exist('S_burst','var') == 1
     ADCP.spec.SSpec.burst=S_burst;
     ADCP.spec.SSpec.units=units;
     ADCP.spec.SSpec.freq=P_freq';
     else
     end
     if exist('V_burst','var') == 1
     ADCP.spec.VSpec.burst=V_burst;
     ADCP.spec.VSpec.units=units;
     ADCP.spec.VSpec.freq=V_freq';
     else
     end
     if exist('D_burst','var') == 1
     ADCP.spec.DSpec.burst=D_burst;
     ADCP.spec.DSpec.units=dir_units;
     ADCP.spec.DSpec.dir_bins=dbinedges;
     ADCP.spec.DSpec.freq=P_freq';
     else
     end
       end
end
ADCP.spec.Time=TIME';
aa = find(ADCP.spec.Time>=stime & ADCP.spec.Time <=etime);
ADCP.spec.time = ADCP.spec.time(aa);
ADCP.spec.PSpec.burst = ADCP.spec.PSpec.burst(:,aa);
ADCP.spec.SSpec.burst = ADCP.spec.SSpec.burst(:,aa);
ADCP.spec.VSpec.burst = ADCP.spec.VSpec.burst(:,aa);
ADCP.spec.DSpec.burst = ADCP.spec.DSpec.burst(:,:,aa);

%========================================================================
 %%
 %====================== LOG 9 TXT file Processing =======================
 % The waves data stored within the log9 file is read in, and th data is
 % stored in the ADCP structure under the field log9
   logfile_source=fullfile(inpath,'*LOG9.TXT');
    logfile_dir=dir(logfile_source);
    file=append(inpath, logfile_dir.name);
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
        grid on; ylabel('D_p (deg)'); axis tight; datetick('x','mm/dd','keeplimits','keepticks')
        
%% Saving Wave and Spectrum data to ADCP structure
% For spectra saving, moved above
 % Data Structure for Log 9 File
 %Trim the beginning and ends of the files:
wave_time2=convertTo(wave_time, 'datenum');
aa = find(wave_time2>=stime & wave_time2 <=etime);
     ADCP.Hs=H_s(aa);
     ADCP.Tp=T_p(aa);
     ADCP.Dp=D_p(aa);
     ADCP.Tp_Sea = Tp_Sea(aa);
     ADCP.Dp_Sea = Dp_Sea(aa);
     ADCP.Hs_Sea = Hs_Sea(aa);
     ADCP.Tp_Swell = Tp_Swell(aa);
     ADCP.Dp_Swell = Dp_Swell(aa);
     ADCP.Hs_Swell = Hs_Swell(aa);
     ADCP.Hmax = Hmax(aa);
     ADCP.Tmax = Tmax(aa);
     ADCP.Honethird = Honethird(aa);
     ADCP.Tonethird = Tonethird(aa);
     ADCP.Honetenth = Honetenth(aa);
     ADCP.Tonetenth = Tonetenth(aa);
     ADCP.Hmean = Hmean(aa);
     ADCP.Tmean = Tmean(aa);
     ADCP.D_pmean=Dmean(aa);
     ADCP.units={'Hs= m'; 'Tp = s'; 'Md = deg'};
     ADCP.Time= wave_time2(aa);



            
        
        
        
               
               