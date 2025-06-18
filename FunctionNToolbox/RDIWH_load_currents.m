function [ADCPc] = RDIWH_load_currents(file_path,save_to_path,depnum,stationID,ADCPparams,stime,etime)  
% [ADCPc] = RDIWH_load_currents(file_path,save_to_path,depnum,stationID,ADCPparams,stime,etime) 
% 
% ======================================================================
   % This function inputs the path of a deployment folder from Wavesmon
   % processing. 
   %
   % CHANGED TO INPUT THE FULL PATH TO THE .PD0 file of interest
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
  wave_source=fullfile(file_path);
  %wave_source=fullfile(file_path,'*.PD0');
  pd0file=dir(wave_source);
  pd0filedir=fullfile(pd0file.folder, pd0file.name);

  % Surface range and ADCP bins below the surface to delete:
surfrange       = ADCPparams.surfrange;
bins2del        = ADCPparams.bins2del;
beamnum         = ADCPparams.beamnum;
qctouse         = ADCPparams.qctouse;
mounting_dist   = ADCPparams.mounting_dist; % Need to check if we actually need this
  
%% 1. Run RDRADCP
  [ADCP] = rdradcp(pd0filedir,1,-1)
  disp('Section 1 complete')

%% Section 2, find the bad values at the end of deployments
% bad values will have an mtime = 0
    index = find(ADCP.mtime == 0);
    ADCP.mtime(index) = [];
    ADCP.number(index)=[];
    ADCP.pitch(index)=[];
    ADCP.roll(index)=[];
    ADCP.heading(index)=[];
    ADCP.pitch_std(index)=[];
    ADCP.roll_std(index)=[];
    ADCP.heading_std(index)=[];
    ADCP.depth(index)=[];
    ADCP.Temperature(index)=[];
    ADCP.salinity(index)=[];
    ADCP.Pressure(index)=[];
    ADCP.Pressure_std(index)=[];
    ADCP.east_vel(:,index)=[];
    ADCP.north_vel(:,index)=[];
    ADCP.vert_vel(:,index)=[];
    ADCP.error_vel(:,index)=[];
    ADCP.corr(:,:,index)=[];
    ADCP.status(:,:,index)=[];
    ADCP.intens(:,:,index)=[];
    ADCP.bt_range(:,index)=[];
    ADCP.bt_vel(:,index)=[];
    ADCP.bt_corr(:,index)=[];
    ADCP.bt_ampl(:,index)=[];
    ADCP.bt_perc_good(:,index)=[];
  disp('Section 2 complete')

%% Section 3 Calculating depth from pressure data
%the provided depth values are discrete so here we calculate our own depth
%field and call it "new_depth" - this will be the preferred depth estimate.

qcnan_dep = NaN(size(ADCP.east_vel));

    if isempty(ADCP.Pressure) == 0
        for ii = 1:length(ADCP.mtime)
            new_depth(ii) = nanmean((ADCP.Pressure(ii))*10/(1026*9.81)); %recorded pressure is in deca-pascals so we use density of water to be 1026 and g to be 9.81
            surf_i = find(abs(ADCP.config.ranges) <= new_depth(ii));
            if isempty(surf_i) == 0 
            surfbin_p(ii) = surf_i(end);
            qcnan_dep(1:surfbin_p(ii)-bins2del,ii) = 0;
            else
            surfbin_p(ii) = NaN;
            end
        end
    else
        new_depth = []; %if there was no valid pressure data the field is empty
        surfbin_p = [];
    end
disp('Section 3 complete')
% Convert those Decapascal units to Dbar:
ADCP.Pressure = ADCP.Pressure./1000;
%% Section 4 Calculating depth from intensity data, called int_depth
qcnan_int = NaN(size(ADCP.east_vel));

for ii = 1:length(ADCP.mtime)
surf_i = find(ADCP.intens(surfrange,beamnum,ii)==max(ADCP.intens(surfrange,beamnum,ii)));
surfbin_i(ii) = surfrange(surf_i(end));
qcnan_int(1:surfbin_i(ii)-bins2del,ii) = 0;
int_depth(ii) = abs(ADCP.config.ranges(surfbin_i(ii)));
end
disp('Section 4 complete')
%% Section 5 Use one of these QC flags to then mask velocity:
if qctouse == 'pres'
    qcflag = qcnan_dep;
    ADCP.depth = new_depth;
    ADCP.qcFlag = ~isnan(qcflag);
elseif qctouse == 'intn'
    qcflag = qcnan_int;
    ADCP.depth = int_depth;
    ADCP.qcFlag = ~isnan(qcflag);
else 
    qcflag = ones(size(ADCP.east_vel));
end

    ADCP.Velocity_East = ADCP.east_vel + qcflag;
    ADCP.Velocity_North = ADCP.north_vel + qcflag;
    ADCP.Velocity_Up = ADCP.vert_vel + qcflag;
    ADCP.bins = abs(ADCP.config.ranges)+mounting_dist;
    ADCP.SN = ADCP.config.remus_serialnum;
    ADCP.Time = ADCP.mtime;
    disp('Section 5 complete')
    
%% Section 6 Remove fields that are not needed for L0
fields = {'name','config','pitch_std','roll_std','heading_std','status','bt_range','bt_vel',...
    'bt_corr','bt_ampl','bt_perc_good','perc_good', 'number', 'Pressure_std','salinity', ...
    'east_vel','north_vel','vert_vel','coords', 'rerunparameter','mtime'};

ADCP = rmfield(ADCP,fields);
disp('Section 6 complete')

ADCP.units={'pres = dB'; 'intens/amp/correl = counts'; 'vels = m/s'};
ADCP.notes = {'Time base is GMT';
    sprintf('Created on %s',datestr(now))};
%%
% Trim the beginning and ends of the files:
aa = find(ADCP.Time>=stime & ADCP.Time <=etime);
ADCPc.Time = ADCP.Time(aa);
ADCPc.pitch = ADCP.pitch(aa);ADCPc.roll = ADCP.roll(aa);ADCPc.heading = ADCP.heading(aa);
ADCPc.Temperature = ADCP.Temperature(aa);ADCPc.Pressure = ADCP.Pressure(aa);ADCPc.depth = ADCP.depth(aa);
ADCPc.bins = ADCP.bins;
ADCPc.units = ADCP.units;
ADCPc.SN = ADCP.SN;
ADCPc.notes = ADCP.notes;

ADCPc.Velocity_East = ADCP.Velocity_East(:,aa);
ADCPc.Velocity_North = ADCP.Velocity_North(:,aa);
ADCPc.Velocity_Up = ADCP.Velocity_Up(:,aa);
ADCPc.velocity_error = ADCP.error_vel(:,aa);
ADCPc.correlation = ADCP.corr(:,:,aa);
ADCPc.amplitude = ADCP.intens(:,:,aa);
ADCPc.qcFlag = ADCP.qcFlag(:,aa);

%% Need the raw data plots here:

%% Saving the file name
%filename=[save_to_path,'SBE_',sprintf('%08d',ADCP.SN),'_','DEP',num2str(depnum),'_',stationID,'_L0.mat'];
%save(filename,'SN','time','pitch','roll','heading','Temperature','Pressure','depth','bins','amplitude','correlation',...
%    'Velocity_East','Velocity_North','velocity_up','velocity_error','qcFlag','units','notes'); 