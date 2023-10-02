function [ADCP, RBRtri, SBE37, RBRsoloT]=L1_BOEMprocessing(deploypath, deploydepth, surfacebin, startTime, endTime)
% L1 BOEM Processing Function
% ========================================================================
% L1_Boemprocessing inputs the file path of the files from the load
% functions to interpolate the RBRtri and SBE37 to be collecting in water data and
% on 5 minute intervals (Sample rate of RDIWH currents). RDIWH Spec and
% log9 files are edited to only display in-water data as well. Inputs the
% to the function are:
% Input:
%   - deploypath: path created by deployname and crateNum to the specific
%   BOEMTest\Crate# folder
%
%   - deploydepth: deployment depth of the lander is necessary for removing data 
% while the lander is physically out of the water and being lowed down to the sea surface. 
% The variable deploydepth in meters is defined by the depth at which the lander 
% is deployed minus a meter to account for tidal fluctuations. 
%
%   - surfacebin: %It defines the sea surface bin # where the echo intensities 
% are > 90 associated with the water depth and bin size. The bin number can be determined
% by looking at the raw output from the RDIWH_load.m function in ADCP.cur.intens variable or 
% divide the mean depth during the deployment by the bin size to get a bin
% number at the surface approximately
% The numeric bin # of this variable is important to set a parameter for NaN data above the sea surface.  
% 
% Output: 
%  - ADCP: structure for RDIWH data (Spectrum, Log9, and Current Data)
%  - RBRtri: structure for RBR tridente data
%  - SBE37: structure for SeaBird 37 data
%  - RBRsoloT: structure for temperature recordings from the RBR Solo T
%
% =======================================================================

%Creating L1 Folder within the specific deployname and Cratenum folders
L1folder=mkdir(deploypath, 'L1processing')

% find the folder with load function .mat files in L0 file
  L0Folder=dir([deploypath, filesep, 'L0*']);
% add folder to inpath
  inpath=[deploypath, filesep, L0Folder(1).name, filesep];


%Loading in the L0.mat file from the saved folder
source=fullfile(inpath, '*.mat');
fileattr=dir(source);

for ii=1:length(fileattr)
    filename=fileattr(ii).name;
    filepath=append(inpath, '\',filename);
    load(filepath);
end

%% Interpreting the RBRtri, SBE37, ADCP currents and wave characteritics to a single time vector
%Converting variables to datetime to visibly see where the different
%instruments and logs start and end with in-water data.
sbTime=datetime(SBE37.time, 'ConvertFrom', 'datenum');
rbrTime=datetime(RBRtri.time, 'ConvertFrom', 'datenum');
specTime=datetime(ADCP.spec.time, 'ConvertFrom', 'datenum');
curTime=datetime(ADCP.cur.mtime, 'ConvertFrom','datenum');
rbrtemp=datetime(RBRsoloT.time, 'ConvertFrom', 'datenum');

%Creating a Start and End time for each instrumnet
    %SeaBird37 (1 min INT)
    SBE37_ind=find(SBE37.time>= datenum(startTime) & SBE37.time<= datenum(endTime));
      SBE37.conductivity = SBE37.conductivity(SBE37_ind);
      SBE37.pressure = SBE37.pressure(SBE37_ind);
      SBE37.temperature = SBE37.temperature(SBE37_ind);
      SBE37.salinity = SBE37.salinity(SBE37_ind);
      SBE37.time =SBE37.time(SBE37_ind);
%       %Interp SBE37 to ADCP currents
%       SBE37.conductivity = interp1(SBE37.time, SBE37.conductivity, ADCP.cur. mtime, 'linear');
%       SBE37.pressure = interp1(SBE37.time, SBE37.pressure, ADCP.cur.mtime, 'linear');
%       SBE37.temperature = interp1(SBE37.time, SBE37.temperature, ADCP.cur.mtime, 'linear');
%       SBE37.salinity = interp1(SBE37.time, SBE37.salinity, ADCP.cur.mtime, 'linear');
%       SBE37.time = interp1(SBE37.time, SBE37.time, ADCP.cur.mtime, 'linear');

  %RBRtri (1 min INT)
    RBRtri_ind=find(RBRtri.time >= datenum(startTime) & RBRtri.time <= datenum(endTime));
      RBRtri.time=RBRtri.time(RBRtri_ind);
      RBRtri.chlorophyll_a=RBRtri.chlorophyll_a(RBRtri_ind);
      RBRtri.FDOM=RBRtri.FDOM(RBRtri_ind);
      RBRtri.turbidity=RBRtri.turbidity(RBRtri_ind);
%       %Interp RBRtri to ADCP Currents
%       RBRtri.turbidity=interp1(RBRtri.time, RBRtri.turbidity, ADCP.cur.mtime, 'linear');
%       RBRtri.chlorophyll_a=interp1(RBRtri.time, RBRtri.chlorophyll_a, ADCP.cur.mtime, 'linear');
%       RBRtri.FDOM=interp1(RBRtri.time, RBRtri.FDOM, ADCP.cur.mtime, 'linear');
%       RBRtri.time=interp1(RBRtri.time, RBRtri.time, ADCP.cur.mtime, 'linear');

    %RBRsoloT 
    if isempty(RBRsoloT) == 0
        RBRtemp_ind=find(RBRsoloT.time>= datenum(startTime) & RBRsoloT.time <= datenum(endTime));
       RBRsoloT.temp=RBRsoloT.temp(RBRtemp_ind);
       RBRsoloT.time=RBRsoloT.time(RBRtemp_ind);
   %Interp RBRsoloT to SBE37 (1 min INT)
       RBRsoloT.temp=interp1(RBRsoloT.time, RBRsoloT.temp, SBE37.time, 'linear');
       RBRsoloT.time=interp1(RBRsoloT.time, RBRsoloT.time, SBE37.time, 'linear');
    end
    
   %ADCP Currents (5 min INT)
   ADCPtemp_ind = find(ADCP.cur.mtime >= datenum(startTime) & ADCP.cur.mtime<= datenum(endTime));
    ADCP.cur.mtime = ADCP.cur.mtime(ADCPtemp_ind);
    ADCP.cur.pitch = ADCP.cur.pitch(ADCPtemp_ind);
    ADCP.cur.roll = ADCP.cur.roll(ADCPtemp_ind);
    ADCP.cur.heading = ADCP.cur.heading(ADCPtemp_ind);
    ADCP.cur.depth = ADCP.cur.depth(ADCPtemp_ind);
    ADCP.cur.temperature = ADCP.cur.temperature(ADCPtemp_ind);
    ADCP.cur.pressure = ADCP.cur.pressure(ADCPtemp_ind);
    ADCP.cur.mean_depth = ADCP.cur.mean_depth(ADCPtemp_ind);
    ADCP.cur.east_vel = ADCP.cur.east_vel(:, ADCPtemp_ind);
    ADCP.cur.north_vel = ADCP.cur.north_vel(:, ADCPtemp_ind);
    ADCP.cur.vert_vel = ADCP.cur.vert_vel(:, ADCPtemp_ind);
    ADCP.cur.error_vel = ADCP.cur.error_vel(:, ADCPtemp_ind);
    ADCP.cur.corr = ADCP.cur.corr(:, :, ADCPtemp_ind);
    ADCP.cur.intens = ADCP.cur.intens(:, :, ADCPtemp_ind); 
    clearvars ADCPtemp_ind

  %ADCP Log9 File (20 min Ensembles- 40 min INT)
  ADCPtemp_ind = find(ADCP.log9.wave_time >= datenum(startTime) & ADCP.log9.wave_time <= datenum(endTime));
    ADCP.log9.Hs = ADCP.log9.Hs(ADCPtemp_ind);
    ADCP.log9.Tp = ADCP.log9.Tp(ADCPtemp_ind);
    ADCP.log9.Dp = ADCP.log9.Dp(ADCPtemp_ind);
    ADCP.log9.Hs_Sea = ADCP.log9.Hs_Sea(ADCPtemp_ind);
    ADCP.log9.Dp_Sea = ADCP.log9.Dp_Sea(ADCPtemp_ind);
    ADCP.log9.Tp_Sea = ADCP.log9.Tp_Sea(ADCPtemp_ind);
    ADCP.log9.Hs_swell = ADCP.log9.Hs_Swell(ADCPtemp_ind);
    ADCP.log9.Tp_Swell = ADCP.log9.Tp_Swell(ADCPtemp_ind);
    ADCP.log9.Dp_Swell = ADCP.log9.Dp_Swell(ADCPtemp_ind);
    ADCP.log9.Hmax = ADCP.log9.Hmax(ADCPtemp_ind);
    ADCP.log9.Tmax = ADCP.log9.Tmax(ADCPtemp_ind);
    ADCP.log9.Honethird = ADCP.log9.Honethird(ADCPtemp_ind);
    ADCP.log9.Tonethird = ADCP.log9.Tonethird(ADCPtemp_ind);
    ADCP.log9.Honetenth = ADCP.log9.Honetenth(ADCPtemp_ind);
    ADCP.log9.Tonetenth = ADCP.log9.Tonetenth(ADCPtemp_ind);
    ADCP.log9.Hmean = ADCP.log9.Hmean(ADCPtemp_ind);
    ADCP.log9.Tmean = ADCP.log9.Tmean(ADCPtemp_ind);
    ADCP.log9.D_pmean = ADCP.log9.D_pmean(ADCPtemp_ind);
    ADCP.log9.wave_time = ADCP.log9.wave_time(ADCPtemp_ind);
    clearvars ADCPtemp_ind;

    %ADCP Spectrum Files (20 min Ensembles- 40 min INT)
    ADCPtemp_ind = find(ADCP.spec.time >= datenum(startTime) & ADCP.spec.time <= datenum(endTime))
    ADCP.spec.DSpec.burst = ADCP.spec.DSpec.burst(:,:,ADCPtemp_ind);
    ADCP.spec.SSpec.burst = ADCP.spec.SSpec.burst(:, ADCPtemp_ind);
    ADCP.spec.VSpec.burst = ADCP.spec.VSpec.burst(:, ADCPtemp_ind);
    ADCP.spec.PSpec.burst = ADCP.spec.PSpec.burst(:, ADCPtemp_ind);
    ADCP.spec.time = ADCP.spec.time(ADCPtemp_ind);


%Determining in and out of water time by the RDIWH recorded depth
for i = 1:length(ADCP.cur)
       if length(unique(ADCP.cur(i).depth))>2 && mean(ADCP.cur(i).depth) > deploydepth %Checking to make sure pressure sensor is recording and in the range of what we expect
           
           bad_depth = find(ADCP.cur(i).depth < 10); %used to index values that are assumed to be above the surface. All values where depth is less than 5 are deleted if they meet the if statement above

           ADCP.cur(i).mtime(bad_depth)=[];
           ADCP.cur(i).pitch(bad_depth)=[];
           ADCP.cur(i).roll(bad_depth)=[];
           ADCP.cur(i).heading(bad_depth)=[];
           ADCP.cur(i).depth(bad_depth)=[];
           ADCP.cur(i).temperature(bad_depth)=[];
           ADCP.cur(i).pressure(bad_depth)=[];
           ADCP.cur(i).east_vel(:,bad_depth)=[];
           ADCP.cur(i).north_vel(:,bad_depth)=[];
           ADCP.cur(i).vert_vel(:,bad_depth)=[];
           ADCP.cur(i).error_vel(:,bad_depth)=[];
           ADCP.cur(i).intens(:,:,bad_depth)=[];
           ADCP.cur(i).mean_depth(bad_depth)=[];
           ADCP.cur(i).corr(:,:, bad_depth) = [];
       end
end

%The echo amplitude is now used to determine pre and post deployment data.
%The echo amplitude (intensity) is expected to be lower in air than in
%water. With a 0.5 meter bin size, Bin 25 is near the surface and has a high intensity, if it is low we
%assume that the instrument is out of the water.
for i = 1:length(ADCP.cur)
   clear bad_echo
   %intens(#, , ) needs to be set at the surface intensities where it begin to spike from the surface
   bad_echo = find(ADCP.cur(i).intens(surfacebin,1,[1:50 end-50:end]) < 90);

   ADCP.cur(i).mtime(bad_echo)=[];
   ADCP.cur(i).pitch(bad_echo)=[];
   ADCP.cur(i).roll(bad_echo)=[];
   ADCP.cur(i).heading(bad_echo)=[];
   ADCP.cur(i).depth(bad_echo)=[];
   ADCP.cur(i).temperature(bad_echo)=[];
   ADCP.cur(i).pressure(bad_echo)=[];
   ADCP.cur(i).east_vel(:,bad_echo)=[];
   ADCP.cur(i).north_vel(:,bad_echo)=[];
   ADCP.cur(i).vert_vel(:,bad_echo)=[];
   ADCP.cur(i).error_vel(:,bad_echo)=[];
   ADCP.cur(i).intens(:,:,bad_echo)=[];
   ADCP.cur(i).mean_depth(bad_echo)=[];
   ADCP.cur(i).corr(: ,:, bad_echo) = [];

%Second form of QC on data by clearing any pressure data less than 2 dB.
% for gg=1:length(SBE37)
%    sbepres=find(SBE37.pressure<= 2);
%    SBE37.pressure(sbepres) = [];
%    SBE37.conductivity(sbepres) = [];
%    SBE37.temperature(sbepres) = [];
%    SBE37.salinity(sbepres) = [];
%    SBE37.time(sbepres) = [];
%   % RBRtri
%    RBRtri.chlorophyll_a(sbepres) = [];
%    RBRtri.FDOM(sbepres) = [];
%    RBRtri.turbidity(sbepres) = [];
%    RBRtri.time(sbepres) = [];
%   % RBRsoloT
%    RBRsoloT.temp(sbepres) = [];
%    RBRsoloT.time(sbepres) = [];
%   % RDI 
%    ADCP.cur.mtime(sbepres) = [];
%    ADCP.cur.pitch(sbepres) = [];
%    ADCP.cur.roll(sbepres) = [];
%    ADCP.cur.heading(sbepres) = [];
%    ADCP.cur.depth(sbepres) = [];
%    ADCP.cur.temperature(sbepres) = [];
%    ADCP.cur.salinity(sbepres) = [];
%    ADCP.cur.pressure(sbepres) = [];
%    ADCP.cur.east_vel(:,sbepres) = [];
%    ADCP.cur.north_vel(:,sbepres) = [];
%    ADCP.cur.vert_vel(:,sbepres) = [];
%    ADCP.cur.error_vel(:,sbepres) = [];
%    ADCP.cur.intens(:,:,sbepres) = [];
%    ADCP.cur.mean_depth(sbepres) = [];
%    ADCP.cur.corr(:, :, sbepres) = [];
% end

%% RDI WH QA/QC procedures to RDIWH
% north values
% velocity error test with a threshold of 0.05 cm/s
for zz=1:length(ADCP.cur)
    [row,col] = find(abs(ADCP.cur(zz).error_vel)>=0.05); %indexes bin and ensemble for which the velocity error is above threshold

    for xx = 1:length(row)
           ADCP.cur(zz).north_vel(row(xx),col(xx))=nan; %changes the indexed bin and ensemble to nan
    end

    %correlation magnitude test with threshold of 110
    for xx=1:length(ADCP.cur(zz).mtime)
        [row,col] = find(ADCP.cur(zz).corr(:,:,xx)<=110); %indexes values that fail
         for yy= 1:length(row)
             ADCP.cur(zz).north_vel(row(yy),xx) = nan; %changes values that failed to nan
         end
    end

    % Echo Amplitude test with a threshold of 10
    %the threshold of 10 referes to the difference in intensity (echo
    %amplitude) between adjacent bins 
    for xx=1:length(ADCP.cur(zz).mtime)
        for ll = 1:4
            for yy=20:ADCP.cur(zz).config.n_cells-1  %This test is only applied to bins above the 15th bin
                dif = ADCP.cur(zz).intens(yy,ll,xx)-ADCP.cur(zz).intens(yy+1,ll,xx);
                ind = find(abs(dif)>=10);

                if sum(ind) >= 1
                   ADCP.cur(zz).north_vel(yy:end,xx)= nan;
                else
                end
            end
        end
    end
end
%% East values
% same as the above procedures only now the east velocity values are
% screened
for zz=1:length(ADCP.cur)
    clear col dif ind row xx yy 
    [row,col] = find(abs(ADCP.cur(zz).error_vel)>=0.05);

    for xx = 1:length(row)
           ADCP.cur(zz).east_vel(row(xx),col(xx))=nan; 
    end

    %correlation magnitude test
    for xx=1:length(ADCP.cur(zz).mtime)
        [row,col] = find(ADCP.cur(zz).corr(:,:,xx)<=110);
         for yy= 1:length(row)
             ADCP.cur(zz).east_vel(row(yy),xx) = nan;
         end
    end

    % Echo Amplitude test
     for xx=1:length(ADCP.cur(zz).mtime)
        for ll = 1:4
            for yy=20:ADCP.cur(zz).config.n_cells-1
                dif = ADCP.cur(zz).intens(yy,ll,xx)-ADCP.cur(zz).intens(yy+1,ll,xx);
                ind = find(abs(dif)>=10);

                if sum(ind) >= 1
                   ADCP.cur(zz).east_vel(yy:end,xx)= nan;
                else
                end
            end
        end
    end
end

%% Vertical Velocity QC
% same as the above procedures only now the vertical velocity values are
% screened
for zz=1:length(ADCP.cur)
    clear col dif ind row xx yy 
    [row,col] = find(abs(ADCP.cur(zz).error_vel)>=0.05);

    for xx = 1:length(row)
           ADCP.cur(zz).vert_vel(row(xx),col(xx))=nan; 
    end

    %correlation magnitude test
    for xx=1:length(ADCP.cur(zz).mtime)
        [row,col] = find(ADCP.cur(zz).corr(:,:,xx)<=110);
         for yy= 1:length(row)
             ADCP.cur(zz).vert_vel(row(yy),xx) = nan;
         end
    end

    % Echo Amplitude test
     for xx=1:length(ADCP.cur(zz).mtime)
        for ll = 1:4
            for yy=20:ADCP.cur(zz).config.n_cells-1
                dif = ADCP.cur(zz).intens(yy,ll,xx)-ADCP.cur(zz).intens(yy+1,ll,xx);
                ind = find(abs(dif)>=10);

                if sum(ind) >= 1
                   ADCP.cur(zz).vert_vel(yy:end,xx)= nan;
                else
                end
            end
        end
    end
end

%% Calculating 4 beam average Echo Amplitude
%Here we average the beam intensity (echo amplitude)
%then we determine a conservative estimate for the sea surface using the
%averaged echo amplitude
if isempty(ADCP.cur.pressure)==1
    disp('pres')
elseif isempty(ADCP.cur.pressure)== 0
for gg = 1:length(ADCP.cur)
    for  hh = 1:length(ADCP.cur(gg).mtime)
         ADCP.cur(gg).EA_avg(:,hh) = mean(ADCP.cur(gg).intens(:,:,hh),2); %averaging across the 4 beams
    end
end

%here we use the same procedure as the echo amplitude QC test. A difference
%of 10 is found between two adjacent bins, except this time the 4 beam
%average echo intensity is used
% for zz = 1:length(ADCP.cur)
%     % create a set of row indices
%     inds = repmat([1:ADCP.cur(zz).config.n_cells]',1,length(ADCP.cur(zz).mtime));
%     % get first difference in rows
%     dum0 = ADCP.cur(zz).EA_avg(1:ADCP.cur(zz).config.n_cells-1,:);
%     dum1 = ADCP.cur(zz).EA_avg(2:ADCP.cur(zz).config.n_cells  ,:);
%     dif  = dum1-dum0;
%     % set the indices of any row with abs(dif)<20 to infinity
%     inds(abs(dif)<20)=inf;
%     c    = abs(min(inds,[],1));

    for xx=1:length(ADCP.cur(zz).mtime)
        for yy=15:ADCP.cur(zz).config.n_cells-1 %only bins above the 15th bin
            dif(xx,yy) = ADCP.cur(zz).EA_avg(yy,xx)-ADCP.cur(zz).EA_avg(yy+1,xx);

            if isempty(min(find(abs(dif(xx,:))>10))) == 0
                c(xx) = min(find(abs(dif(xx,:))>10));
            else
                c(xx) = 25;
            end
        end
    end

%the config ranges is used to go from bin number to distance from ADCP
%(depth). then a mean across the entire deployment is taken
   ADCP.cur(zz).depth_ea = mean(abs(ADCP.cur(zz).config.ranges(c))); 
end

%% Removing values that are above Sea surface 
%last resort QC procedure
%anything above the sea surface as determined by depth is turned into an
%nan
%%%%
%if there is no depth field the sea surface determined by the echo
%amplitude is used as a conservative cut off

for ii=1:length(ADCP.cur)
    if isempty(ADCP.cur.mean_depth) == 0

        surface_press = find(abs(ADCP.cur.config.ranges) > ADCP.cur.mean_depth);
        surface_cut_press = min(surface_press);

        ADCP.cur.east_vel(surface_cut_press:end,:) = nan;
        ADCP.cur.north_vel(surface_cut_press:end,:) = nan;
        ADCP.cur.vert_vel(surface_cut_press:end,:) = nan;

    else
        surface_ea = find(abs(ADCP.cur.config.ranges) > ADCP.cur.depth_ea);
        surface_cut_ea = min(surface_ea);

        ADCP.cur.east_vel(surface_cut_ea:end,:) = nan;
        ADCP.cur.north_vel(surface_cut_ea:end,:) = nan;
        ADCP.cur.vert_vel(surface_cut_ea:end,:) = nan;
    end
end
disp('QA/QC Done')

%% Saving File to L1 processing 
date=datestr(ADCP.spec.time(1),'yyyy_mm_dd');
L1folder=dir([deploypath, filesep, 'L1pr*']);
%Changing dir back to BOEMTest folder
outpath=[inpath, '..',]; cd(outpath); outpath=pwd;
filename=[outpath, filesep,L1folder(1).name, filesep, 'BOEML1proc_', date,'.mat'];
save(filename, 'ADCP','RBRtri','SBE37', 'RBRsoloT'); 
    
end