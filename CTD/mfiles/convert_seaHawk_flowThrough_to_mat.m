function seahawk = convert_seaHawk_flowThrough_to_mat(datDir,dates,arcDir);
%
% USAGE: seahawk = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);
%
% inputs:
% datDir= top directory for raw data files
% date  = {'yyyymmdd'} cell array of sub-directories 
% arcDir= directory to archive data (Default: datDir)
% outputs:
% seahawkFT = structure array containing flow-through data.

%
if nargin == 3
    outDir = arcDir;
else
    outDir = datDir;
end
%
nd        = length(dates);
temporary = struct([]);
seahawk   = struct([]);
NMEA      = struct([]);
SENSOR    = struct([]);
for jj = 1:nd
    date = dates{jj};
    % Data directory and data files structure
    flowDir   = [datDir,date];
    flowFiles = dir(flowDir);
    % remove '.', '..', and '.mat' files
    flowFiles(ismember({flowFiles.name},{'.','..'})) = [];
    flowFiles(endsWith({flowFiles.name},'.mat')) = [];
    allFiles  = {flowFiles.name}';
    % find indices of NMEA & flow-through files 
    NMEAexp   = regexp(allFiles,'.*(NMEA).*');
    isNMEA    = ~cellfun('isempty',NMEAexp);
    RVSHexp   = regexp(allFiles,'.*(RVSH).*');
    isRVSH    = ~cellfun('isempty',RVSHexp);
    % split files
    if sum(isNMEA)~=sum(isRVSH), disp('missing at least one NMEA/RVSH file pair'), return, end
    NMEAfiles = flowFiles(isNMEA);
    RVSHfiles = flowFiles(isRVSH);    
    nf = sum(isNMEA);
    % process the NMEA data first and match RVSH to closest day/time NMEA
    for kk = 1:nf
        NMEAfile = [NMEAfiles(kk).folder, filesep, NMEAfiles(kk).name];
        fid      = fopen(NMEAfile);
        formats  = {'GPGGA', '%*s %f %f %s %f %s %f %f %f %f %s %f %s %f %*s';...
                    'GPGLL', '%*s %f %s %f %s %f %s %s %*s';...
                    'GPVTG', '%*s %f %s %f %s %f %s %f %s %*s';...
                    'GLZDA', '%*s %f %d %d %d %s %*s'};
        % start w/ a GPGGA string
        nLines   = 0;
        flag     = 0;
        while ~feof(fid)% Fnum>1
            tline = fgetl(fid);
            nLines = nLines+1;
            Fnum  = find(ismember(formats(:,1),tline(2:6)));
            if Fnum==1 & ~flag
                startPointer = ftell(fid)-(length(tline)+2);
                flag         = 1;
            end
        end
        % rewind to start of current line
        fseek(fid,startPointer,-1);
        % initialize variables
        hour0 = nan(nLines,1); min0 = hour0; sec0 = hour0;
        hour1 = nan(nLines,1); min1 = hour1; sec1 = hour1;        
        year  = nan(nLines,1); mon  = year ; day  = year ;
        lat0  = nan(nLines,1); lon0 = lat0 ; lat1 = lat0 ; lon1 = lat0;
        % zero iterating vars
        Fnum     = 0;
        lnum     = 0;
        while ~feof(fid)
% $$$             if lnum==3; break; end            
            tline = fgetl(fid);
            %
            Fnew = find(ismember(formats(:,1),tline(2:6)));
            if isempty(Fnew)
                continue
            end
% $$$             if (isempty(Fnew) | (Fnew~=Fnum+1)) %& ~flag
% $$$                 Fnum=0;
% $$$                 %                flag=1;
% $$$                 continue
% $$$             elseif isempty(Fnew) %| flag
% $$$                 continue
% $$$             else
% $$$                 flag=0;
% $$$             end
            %
            % parse string
            D = textscan(tline,formats{Fnew,2},'Delimiter',',');
            switch Fnew
              case 1
                lnum = lnum+1;
                % GPGGA
                % time hhmmss.ss
                hhmmss = D{1};
                ss     = round(rem(hhmmss,1e2),4);
                mm     = rem(hhmmss-ss,1e4)/1e2;
                hh     = (hhmmss-mm*1e2-ss)/1e4;
                hour0(lnum) = hh;
                min0(lnum)  = mm;
                sec0(lnum)  = ss;
                % lat ddmm.mmmm
                ddmm   = D{2};
                NS     = (char(D{3})=='N')-0.5;
                mm     = round(rem(ddmm,100),6);
                lat0(lnum) = sign(NS)*((ddmm-mm)/100 + mm/60);
                % lon ddmm.mmmm
                ddmm   = D{4};
                EW     = (char(D{5})=='E')-0.5;
                mm     = round(rem(ddmm,100),6);
                lon0(lnum) = sign(EW)*((ddmm-mm)/100 + mm/60);
                fix(lnum) = D{6};
                nSat(lnum)= D{7};
                hdop(lnum)= D{8};
              case 2

                % lat ddmm.mmmm
                ddmm   = D{1};
                NS     = (char(D{2})=='N')-0.5;
                mm     = round(rem(ddmm,100),6);
                lat1(lnum) = sign(NS)*((ddmm-mm)/100 + mm/60);
                % lon ddmm.mmmm
                ddmm   = D{3};
                EW     = (char(D{4})=='E')-0.5;
                mm     = round(rem(ddmm,100),6);
                lon1(lnum) = sign(EW)*((ddmm-mm)/100 + mm/60);
                status(lnum) = D{5};
                % time hhmmss.ss
                hhmmss = D{1};
                ss     = round(rem(hhmmss,1e2),4);
                mm     = rem(hhmmss-ss,1e4)/1e2;
                hh     = (hhmmss-mm*1e2-ss)/1e4;
                hour1(lnum) = hh;
                min1(lnum)  = mm;
                sec1(lnum)  = ss;
              case 3
                % GPVTG
                heading(lnum) = D{1};
                headRef(lnum) = char(D{2});
                speed(lnum)   = D{7}*1e3/3600;% km/hr --> m/s
                speedRef(lnum,:)= 'mps';
              case 4
                % GPZDA utc time
                % time hhmmss.ss
                hhmmss = D{1};
                ss     = round(rem(hhmmss,1e2),4);
                mm     = rem(hhmmss-ss,1e4)/1e2;
                hh     = (hhmmss-mm*1e2-ss)/1e4;
% $$$                 hour(lnum) = hh;
% $$$                 min(lnum)  = mm;
% $$$                 sec(lnum)  = ss;
                day(lnum)  = D{2};
                mon(lnum)  = D{3};
                year(lnum) = D{4};
                if ~isempty(D{5}) & ~exist('zone','var')
                   zone    = char(D{5});
                end
            end
% $$$             Fnum = Fnew;
% $$$             if Fnum==4;
% $$$ % $$$                 lnum=lnum+1;
% $$$                 Fnum=0;
% $$$             end
        end
        % trim pre-allocated vars
        N = min(length(speed),lnum);
        year = year(1:N); mon = mon(1:N); day = day(1:N);
        hour0=hour0(1:N); min0=min0(1:N); sec0=sec0(1:N);        
        lat0 =lat0(1:N) ; lon0= lon0(1:N);
        %
        %
        NMEA(jj,kk).date_vec  = [ year, mon, day, hour0, min0, sec0 ];
        NMEA(jj,kk).timeUTC   = datenum([ year, mon, day, hour0, min0, sec0 ]);
        NMEA(jj,kk).time_zone = zone';
        NMEA(jj,kk).latitude  = lat0;
        NMEA(jj,kk).longitude = lon0;
        NMEA(jj,kk).hhmmss    = [ hour0, min0, sec0 ];
        NMEA(jj,kk).heading   = heading';
        NMEA(jj,kk).headRef   = headRef';
        NMEA(jj,kk).speed     = speed';
        NMEA(jj,kk).speedRef  = speedRef;
        NMEA(jj,kk).HDOP      = hdop';
        NMEA(jj,kk).fix       = double(fix');
        NMEA(jj,kk).satellites= double(nSat');
        NMEAstartTimeUTC(kk)  = min(NMEA(jj,kk).timeUTC);
        NMEAendTimeUTC(kk)    = max(NMEA(jj,kk).timeUTC);        
% $$$         figure, geoscatter(lat0,lon0,4,speed','filled'), geobasemap('satellite')        
    end
    %
    %
    for kk = 1:nf
        RVSHfile = [RVSHfiles(kk).folder, filesep, RVSHfiles(kk).name];
        fid      = fopen(RVSHfile);
        flag     = 0;
        while ~flag
            tline = fgetl(fid);
            dum0  = regexp(tline,'.*(MEAN VALUE).*');
            dum1  = regexp(tline,'.*(STANDARD DEVIATION:).*');            
            dum2  = regexp(tline,'.*(Date).*');
            if dum0
                meanString = tline;
            elseif dum1
                stdString = tline;
            elseif dum2
                varString = tline;
                flag  = 1;
            end
        end
        format = '%s %s %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
        D      = textscan(fid,format,'Delimiter',',');
        N      = size(D{1},1);
        timeStr= cat(2, char(D{1}),repmat(' ',N,1),char(D{2}));
%    looks like the time is already in UTC, don't need to add: +ET2UTC/24
% $$$         monthOfSample = month(datetime(timeStr(1,:)));
% $$$         weekOfMonth   = week(datetime(timeStr(1,:)),'weekofmonth');
% $$$         ET2UTC = 5;
% $$$         if (monthOfSample>4 & monthOfSample<11) | (monthOfSample==4 & weekOfMonth>2) | (monthOfSample==11 & weekOfMonth==1)
% $$$             ET2UTC = 4;
% $$$         end
        SENSOR(jj,kk).timeUTC    = datenum(timeStr)+D{3}/86400;
        SENSOR(jj,kk).temperature= D{23};
        SENSOR(jj,kk).salinity   = D{15};
        SENSOR(jj,kk).turbidity  = D{19};
        SENSOR(jj,kk).ChlA       = D{5};
        SENSOR(jj,kk).fDOM       = D{8};
        SENSOR(jj,kk).O2sat      = D{11};
        SENSOR(jj,kk).O2mg_L     = D{13};
        SENSOR(jj,kk).pH         = D{22};        
        SENSOR(jj,kk).pressure   = D{14};
        SENSOR(jj,kk).fileName   = RVSHfiles(kk).name;
        %
        vars                   = split(varString,',');
        sensorColumns          = [23,15,19,5,8,11,13,22,14];
        if length(vars)<=23
            fprintf('\n missing data columns for: %s\n', RVSHfiles(kk).name)
            sensorColumns(sensorColumns>length(vars))=[];
        end
        SENSOR(jj,kk).units    = vars(sensorColumns);
        SENSORstartTimeUTC(kk) = min(SENSOR(jj,kk).timeUTC);
        SENSORendTimeUTC(kk)   = max(SENSOR(jj,kk).timeUTC);        
    end
    %
    if nf>1
        % now pair up the NMEA/SENSOR structures based on start-times
        pairs     = 1:nf;
        startDiff = abs(NMEAstartTimeUTC'-SENSORstartTimeUTC)*86400;
        endDiff   = abs(NMEAendTimeUTC'-SENSORendTimeUTC)*86400;
        % get the index of the nearest NMEA to each SENSOR
        [minStartDiff,indNMEAstart] = min(startDiff,[],1);
        [minEndDiff  ,indNMEAend  ] = min(endDiff  ,[],1);
        startUncertainty       = minStartDiff./min(startDiff(pairs'~=indNMEAstart));
        endUncertainty         = minEndDiff./min(  endDiff(  pairs'~=indNMEAend  ));
        fprintf('fractional error btwn NMEA/SENSOR start-time: %f \n',round(mean(startUncertainty),5))
        indNMEA = indNMEAstart;
    else
        indNMEA=1;
    end
    %
    % combine each record (need to decide which variables to interpret NMEA or SENSOR)
    for kk = 1:nf
        ii         = indNMEA(kk);
        % get times and acceptable limits
        timeSENSOR = SENSOR(jj,kk).timeUTC;
        timeNMEA   = NMEA(jj,ii).timeUTC;
        timeStart  = max( timeSENSOR(1)  , timeNMEA(1)  );
        timeEnd    = min( timeSENSOR(end), timeNMEA(end));
        % interpolate to nearest 1-second.
        timeINT    = timeNMEA(timeNMEA>=timeStart & timeNMEA<=timeEnd);
        timeINT    = round( timeINT*86400 );
        timeINT    = unique(timeINT);
        %
        temporary(jj,kk).time = timeINT/86400;
        temporary(jj,kk).date_vec = datevec(timeINT/86400);
        %
        varsNMEA   = {'latitude','longitude','speed','HDOP','fix','satellites'};
        goodNMEA   = find(diff(timeNMEA*86400)>0);
        for ll = 1:length(varsNMEA)
            eval(['temporary(jj,kk).',varsNMEA{ll},' = interp1(timeNMEA(goodNMEA)*86400,NMEA(jj,ii).',varsNMEA{ll},'(goodNMEA,:),timeINT);']);
        end
        %
        temporary(jj,kk).sensor_units = SENSOR(jj,kk).units;
        varsSENSOR = {'temperature','salinity','turbidity','ChlA','fDOM','O2sat','O2mg_L','pH','pressure'};
        goodSENSOR   = find(diff(timeSENSOR*86400)>0);
        for ll = 1:length(varsSENSOR)
            eval(['temporary(jj,kk).',varsSENSOR{ll},' = interp1(timeSENSOR(goodSENSOR)*86400,SENSOR(jj,ii).',varsSENSOR{ll},'(goodSENSOR,:),timeINT);']);
        end
        %
        % archive one file for each input file
        seahawk = temporary(jj,kk);
        if isempty(seahawk(1).time)
            fprintf('\n empty data structure for: %s\n', SENSOR(jj,kk).fileName)
            continue
        end
        fileOut = sprintf([arcDir,filesep,'seahawk_flowThrough_%s_L0.mat'],datestr(seahawk(1).time(1),'yyyymmdd_HHMM'));
        save(fileOut,'-struct','seahawk')
    end
end
% now save a combined file
seahawk = temporary;
fileOut = sprintf([arcDir,filesep,'seahawk_flowThrough_%s_to_%s_L0.mat'],datestr(seahawk(1).time(1),'yyyymmdd'),datestr(seahawk(end).time(1),'yyyymmdd'));
save(fileOut,'seahawk')
end
