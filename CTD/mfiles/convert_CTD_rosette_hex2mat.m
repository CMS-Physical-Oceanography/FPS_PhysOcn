function [SBE] = convert_CTD_rosette_hex2mat(datDir,date,arcDir);
%
%
% USAGE: [SBE] = convert_CTD_rosette_hex2mat(datDir,date,arcDir);
%
% inputs:
% datDir= top directory for raw data files
% date  = 'yyyymmdd' sub-directory of sample
% arcDir= directory to archive data (Default: datDir)
% outputs:
% SBE   = structure array containing each cast for this date.


if nargin == 3
    outDir = arcDir;
else
    outDir = datDir;
end
%
nd = length(dates);
SBE= struct([]);
for jj = 1:nd
    date = dates{jj};
    % Data directory and data files structure
    CTDdir   = [datDir,date];
    ctdFiles = dir([CTDdir,filesep,'station*.hex']);
    conFile  = dir([CTDdir,filesep,'*.xmlcon']);
    conFile  = [conFile.folder,filesep,conFile.name];
    config   = parseXML(conFile);
    %
    nf     = size(ctdFiles,1);
    castNum= 0;
    %
    for kk = 1:nf
        clear zgrid Pgrid Dgrid D2grid Sgrid Tgrid t x y
        %
        % READ CTD - SBE25plus data file
        % ttttttccccccspppuuuvvv
        %
        % tttttt = temperature frequency, 3 bytes
        % cccccc = conductivity frequency, 3 bytes
        % s = sign character for pressure (sign character 0 = +; sign character 4 = -)
        % ppp = pressure, 12 bits
        % uuu = stored voltage output 0, 12 bits
        % vvv = stored voltage output 1, 12 bits
        fid = fopen([ctdFiles(kk).folder,filesep,ctdFiles(kk).name]);
        pasthdr = 0;
        if fid ~= -1
            fprintf(['processing file: ',ctdFiles(kk).name,' \n'])
            lnum = 0;
            while ~feof(fid)
                %Grab line from file
                tline = fgetl(fid);
                %Look for Lat Lon and UTC time in file header 
                if lnum==0
                    split = regexp(tline,'=','split');
                    str   = split{1};
                    switch str
                      case '* FileName '
                        SBE(kk).rawFile = split{2};
                      case '* Temperature SN '
                        SBE(kk).TemperatureSN = str2num(split{2});                    
                      case '* Conductivity SN '
                        SBE(kk).ConductivitySN = str2num(split{2});
                      case '* NMEA Latitude '
                        dum = strtrim(split{2});
                        pm  = 2*strcmp(dum(end),'N')-1;
                        SBE(kk).Latitude = sign(pm)*(str2num(dum(1:3)) + str2num(dum(4:end-1))/60);
                      case '* NMEA Longitude '
                        dum = strtrim(split{2});
                        pm  = 2*strcmp(dum(end),'E')-1;
                        SBE(kk).Longitude = sign(pm)*(str2num(dum(1:3)) + str2num(dum(4:end-1))/60);
                      case '* NMEA UTC (Time) '
                        dum = strtrim(split{2});
                        SBE(kk).Time_UTC  = datestr(dum);
                    end
                end
                %
                %Begin data collection at end of header
                if min(tline(1:5) == '*END*')
                    pasthdr = 1;
                    lnum    = 1;
                    continue
                end
                
                %Grab data
                if pasthdr
                    % break
                    temphex = tline(1:8);
                    condhex = tline(9:16);
                    preshex  = tline(17:22);
                    ptmphex  = tline(23:28);               
                    aux0hex  = tline(29:32); 
                    aux1hex  = tline(33:36);
                    aux2hex  = tline(37:40);
                    aux3hex  = tline(41:44);
                    aux4hex  = tline(45:48);
                    latihex  = tline(49:54);                              
                    longhex  = tline(55:60);               
                    timehex  = tline(61:70);
                    t0hex    = timehex(9:10);
                    t1hex    = timehex(7:8);
                    t2hex    = timehex(5:6);
                    t3hex    = timehex(3:4);
                    t4hex    = timehex(1:2);
% $$$                aux5hex  = tline(49:52); % not used
% $$$                aux6hex  = tline(53:56); % not used
% $$$                aux7hex  = tline(57:60); % not used
% $$$                diaghex  = tline(61:68); % not used              
                    % algorithm for SBE25plus--needs verification!
                    q32           = quantizer('float','nearest','saturate',[32 8]);
                    tempbin  =  dec2bin(hex2dec(temphex'),4)';
                    tempbin  = (tempbin(:)');
                    tempraw(lnum) = hex2num(q32,temphex);
                    condraw(lnum) = hex2num(q32,condhex);
                    presraw(lnum) = hex2dec(preshex);
                    ptmpraw(lnum) = hex2dec(ptmphex);
                    aux0raw(lnum) = hex2dec(aux0hex)*5.000/2^16;
                    aux1raw(lnum) = hex2dec(aux1hex)*5.000/2^16;
                    aux2raw(lnum) = hex2dec(aux2hex)*5.000/2^16;
                    aux3raw(lnum) = hex2dec(aux3hex)*5.000/2^16;
                    aux4raw(lnum) = hex2dec(aux4hex)*5.000/2^16;
% $$$                aux5raw(lnum) = hex2dec(aux5hex)*5.000/2^16;% not used
% $$$                aux6raw(lnum) = hex2dec(aux6hex)*5.000/2^16;% not used
% $$$                aux7raw(lnum) = hex2dec(aux7hex)*5.000/2^16;% not used
                    %
                    % Convert NMEA stream (custom, determined by flipping bits)
                    lon(lnum) = 0-hex2dec(longhex)*2e-5;
                    lat(lnum) =   hex2dec(latihex)*2e-5;               
                    t0 = hex2dec(t0hex)*2^24;
                    t1 = hex2dec(t1hex)*2^16;
                    t2 = hex2dec(t2hex)*2^8;
                    t3 = hex2dec(t3hex)*2^0;
                    t4 = hex2dec(t4hex)*2^-8;
                    seconds(lnum) = (t0+t1+t2+t3+t4);
                    t(lnum) = datenum('Jan 1 2000') + seconds(lnum)/86400;
                    %
                    lnum = lnum+1;
                end
            end
        end
         fclose(fid);
         %
         rawData = struct('temp',tempraw,'pres',presraw,'ptmp',ptmpraw,'aux0',aux0raw,'aux1',aux1raw,'aux2',aux2raw,'aux3',aux3raw,'aux4',aux4raw);
         % now do conversions to physical units
         break
         data = convert_units_SBE25(rawData,config);
         SBE(jj,castNum) = data;
         castNum = catNum+1;
    end
    % now save level0 data and make a few plots
    figure, geoscatter(lat,lon)
    geobasemap('satellite')
end
% separate into up/down casts?
% process them separately

% $$$  OLD CODE
% $$$ %% get rid of data when CTD is sitting still or ascending
% $$$ 
% $$$ Pu =flip(P);% for upcast
% $$$ Du =flip(D); 
% $$$ Tu =flip(T); 
% $$$ Su =flip(S); 
% $$$ 
% $$$ for ll = 1:length(P)
% $$$     if any(P(ll)<=P(1:ll-1))
% $$$         P(ll) = NaN;
% $$$         D(ll) = NaN;
% $$$         T(ll) = NaN;
% $$$         S(ll) = NaN;
% $$$ % $$$                 Draw(ll) = nan;
% $$$ % $$$                 d0(ll) = nan;
% $$$ % $$$                 d1(ll) = nan;
% $$$ % $$$                 d2(ll) = nan;
% $$$ % $$$                 d3(ll) = nan;
% $$$ % $$$                 d4(ll) = nan;
% $$$ % $$$                 d5(ll) = nan;
% $$$ % $$$                 d6(ll) = nan;
% $$$     end
% $$$ end
% $$$ 
% $$$ for ll = 1:length(P)
% $$$     if any(Pu(ll)<=Pu(1:ll-1))
% $$$         Pu(ll) = NaN;
% $$$         Du(ll) = NaN;
% $$$         Tu(ll) = NaN;
% $$$         Su(ll) = NaN;
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ %% flag dye spikes. 
% $$$ i0 = find(D>120);% don't expect to see this much dye
% $$$ i1 = find(abs(diff(D))>1.5);% gradient of 12ppb/s? 
% $$$ i2 = find(abs(diff(D,2,2))>.75);% 48ppb/s^2?
% $$$ 
% $$$ i3 = union(i0,i1+1);
% $$$ i4 = union(i3,i2+1);
% $$$ 
% $$$ D(i4) = nan;
% $$$ 
% $$$ %% nan all data with bad P
% $$$ inan = find(isnan(P) | isnan(D));
% $$$ P(inan) = [];
% $$$ D(inan) = [];
% $$$ T(inan) = [];
% $$$ S(inan) = [];
% $$$ 
% $$$ %% nan all data with bad P
% $$$ inan = find(isnan(Pu) | isnan(Du));
% $$$ Pu(inan) = [];
% $$$ Du(inan) = [];
% $$$ Tu(inan) = [];
% $$$ Su(inan) = [];
% $$$ 
% $$$ %% log maximum pressure for P-grid later on
% $$$ if kk==1
% $$$     Pmax = nanmax(P);
% $$$ elseif nanmax(P)>Pmax
% $$$     Pmax = nanmax(P);
% $$$ end
% $$$ 
% $$$ 
% $$$ %% fill cell with profiles
% $$$ Dc{kk} = D;
% $$$ Tc{kk} = T;
% $$$ Pc{kk} = P;
% $$$ Sc{kk} = S;
% $$$ 
% $$$ %% fill cell with profiles
% $$$ Duc{kk} = flip(Du);
% $$$ Tuc{kk} = flip(Tu);
% $$$ Puc{kk} = flip(Pu);
% $$$ Suc{kk} = flip(Su);
% $$$ 
% $$$      else 
% $$$          disp(['could not open file: ', [CTDdir,CASTdir,castfiles(indices(kk)).name]])
% $$$          return
% $$$ end 
% $$$ end   
% $$$ 
% $$$ %% interpolate onto P-grid for nanmeans of T(z),D(z),S(z)
% $$$ Pgrid = [0:0.1:ceil(Pmax)]';
% $$$ 
% $$$ 
% $$$ for k = 1:length(Dc)
% $$$     Dgrid(:,k) = [interp1(Pc{k},Dc{k},Pgrid)]';
% $$$     Tgrid(:,k) = [interp1(Pc{k},Tc{k},Pgrid)]';
% $$$     Sgrid(:,k) = [interp1(Pc{k},Sc{k},Pgrid)]';
% $$$     Dugrid(:,k) = [interp1(Puc{k},Duc{k},Pgrid)]';
% $$$     Tugrid(:,k) = [interp1(Puc{k},Tuc{k},Pgrid)]';
% $$$     Sugrid(:,k) = [interp1(Puc{k},Suc{k},Pgrid)]';
% $$$     
% $$$ end
% $$$ 
% $$$ 
% $$$ %% convert lat, lon, time (UTC) --> IB09 x, y, t (pdt)
% $$$ for k = 1:size(latstr,2)
% $$$     lat(k,1) = str2num(latstr{k}(1:2))+ ...
% $$$         str2num(latstr{k}(4:8))/60; % N
% $$$     lon(k,1) = -(str2num(lonstr{k}(1:3))+ ...
% $$$                  str2num(lonstr{k}(5:9))/60);  % W --> negative
% $$$     
% $$$     time0 = datenum('08/31/2015 00:00:00');
% $$$     year = 2015;
% $$$     day = (rem(sample_date(ii),10000) - 15)./100;
% $$$     month = (sample_date(ii) - day)/10000;
% $$$     hour = str2num(tstr_UTC{k}(1:2)) - 7;% PDT = UTC-7
% $$$     mins = str2num(tstr_UTC{k}(4:5));
% $$$     sec = str2num(tstr_UTC{k}(7:8));
% $$$     time = datenum([2015 month day hour mins sec]);
% $$$     t(k,1)= time-time0;
% $$$     t2(k,1) = t_UTC{k} - time0 - 7/24;
% $$$ % $$$         t_UTC(k,1) = 3600*str2num(tstr_UTC{k}(1:2))+60* ...
% $$$ % $$$                           str2num(tstr_UTC{k}(4:5))+ ...
% $$$ % $$$                           str2num(tstr_UTC{k}(7:8));
% $$$ % $$$         t(k,1) = t_UTC(k) - 7*3600; % PDT = UTC - 7 hours
% $$$ end
% $$$ 
% $$$ %% depth approximation [m]
% $$$ zgrid = -Pgrid;
% $$$ 
% $$$ 
% $$$ %% convert from lat lon to xy in IB09 coords
% $$$ [x,y] = lltoxy_imperialbeach(lat,lon);
% $$$ 
% $$$ 
% $$$ %% Compile data into struct
% $$$ SBE(ii).date = datestr(datenum(sprintf('%06d', ...
% $$$                                        sample_date(ii)),'mmddyy')); 
% $$$ SBE(ii).t = t2;
% $$$ SBE(ii).P = Pgrid;
% $$$ SBE(ii).T = Tgrid;
% $$$ SBE(ii).S = Sgrid;
% $$$ SBE(ii).D = Dgrid;
% $$$ SBE(ii).Tu = Tugrid;
% $$$ SBE(ii).Su = Sugrid;
% $$$ SBE(ii).Du = Dugrid;
% $$$ 
% $$$ SBE(ii).lat=lat;
% $$$ SBE(ii).lon=lon;
% $$$ SBE(ii).x = x;
% $$$ SBE(ii).y = y;
% $$$ SBE(ii).z = zgrid;


end
