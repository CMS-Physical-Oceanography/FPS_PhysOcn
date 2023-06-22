function seahawkFT = convert_seaHawk_flowThrough_to_mat(datDir,dates,arcDir);
%
% USAGE: seahawkFT = convert_CTD_rosette_hex2mat(datDir,dates,arcDir);
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
seahawkFT = struct([]);
for jj = 1:nd
    date = dates{jj};
    % Data directory and data files structure
    flowDir   = [datDir,date];
    flowFiles = dir(flowDir);
    allFiles  = {flowFiles.name}';
    % find indices of NMEA & flow-through files 
    NMEAexp   = regexp(allFiles,'.*(NMEA).*');
    isNMEA    = ~cellfun('isempty',NMEAexp);
    RVSHexp   = regexp(allFiles,'.*(RVSH).*');
    isRVSH    = ~cellfun('isempty',RVSHexp);
    % split files
    if sum(isNMEA)~=sum(isRVSH), disp('missing at least one NMEA/RVSH file pair'), regurn, end
    NMEAfiles = flowFiles(isNMEA);
    RVSHfiles = flowFiles(isRVSH);    
    nf = sum(isNMEA);
    % process the NMEA data first and match RVSH to closest day/time NMEA
    for kk = 1:nf
        NMEAfile = [NMEAfiles(kk).folder, filesep, NMEAfiles(kk).name];
        fid      = fopen(NMEAfile);
        formats  = {'GPGGA', '%*s %f %f %s %f %s %u %u %f %f %s %f %s %f %*s';...
                    'GPGLL', '%*s %f %s %f %s %f %s %s %*s';...
                    'GPVTG', '%*s %f %s %f %s %f %s %f %s %*s';...
                    'GLZDA', '%*s %f %d %d %d %s %*s'};
        % make sure that we start w/ a GPGGA string
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
        hour  = nan(nLines,1); min  = hour ; sec  = hour ;
        lat0  = nan(nLines,1); lon0 = lat0 ; lat1 = lat0 ; lon1 = lat0;
        % zero iterating vars
        Fnum     = 0;
        lnum     = 1;
        while ~feof(fid)
            if lnum==26806; break; end            
            tline = fgetl(fid);
            %
            Fnew = find(ismember(formats(:,1),tline(2:6)));
            if (isempty(Fnew) | (Fnew~=Fnum+1)) & ~flag
                Fnum=0;
                flag=1;
                continue
            elseif isempty(Fnew) & flag
                continue
            else
                flag=0;
            end
            %
            % parse string
            D = textscan(tline,formats{Fnew,2},'Delimiter',',');
            switch Fnew
              case 1
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
                unit(lnum)= D{9};
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
                hour(lnum) = hh;
                min(lnum)  = mm;
                sec(lnum)  = ss;
                day(lnum)  = D{2};
                mon(lnum)  = D{3};
                year(lnum) = D{4};
                if ~isempty(D{5}) & ~exist('zone','var')
                   zone(lnum) = char(D{5});
                end
            end
            Fnum = Fnew;
            if Fnum==4;
                lnum=lnum+1;
                Fnum=0;
            end
        end
        % trim pre-allocated vars
        year = year(1:lnum); mon = mon(1:lnum); day = day(1:lnum);
        hour = hour(1:lnum); min = min(1:lnum); sec = sec(1:lnum);
        hour1=hour1(1:lnum); min1=min1(1:lnum); sec1=sec1(1:lnum);
        lat0 =lat0(1:lnum) ; lon0= lon0(1:lnum);
        lat1 =lat1(1:lnum) ; lon1= lon1(1:lnum);
        figure, geoscatter(lat1,lon1), geobasemap('satellite')        
        %
        seahawkFT(jj,kk).date_vec = [ year, mon, day, hour, min, sec ];
        seahawkFT(jj,kk).time     = datenum([ year, mon, day, hour0, min0, sec0 ]);
        seahawkFT(jj,kk).time_zone= zone';
        seahawkFT(jj,kk).lat      = lat0
        seahawkFT(jj,kk).lon      = lon0
        seahawkFT(jj,kk).hhmmss   = [ hour0, min0, sec0 ];
        seahawkFT(jj,kk).heading  = heading';
        seahawkFT(jj,kk).headRef  = headRef';
        seahawkFT(jj,kk).speed    = speed';
        seahawkFT(jj,kk).speedRef = speedRef;
        seahawkFT(jj,kk).HDOP     = hdop;
        seahawkFT(jj,kk).fix      = fix;
        seahawkFT(jj,kk).num_sat  = nSat;
        seahawkFT(jj,kk).unit     = unit;
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
        format = '%s %s %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
        D      = textscan(fid,format,'Delimiter',',');
        N      = size(D{1},1);
        timeStr= cat(2, char(D{1}),repmat(' ',N,1),char(D{2}));
        timeUTC= datenum(timeStr)+4;
    end
end
end