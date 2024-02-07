files = dir('../data/*.rsk');
% $$$ startStop = [datenum('07-Feb-2024 20:00:00'); datenum('07-Feb-2024 21:50:00')];
% $$$ useStartStop = 1;
for ii=1:length(files)
    %
    fin = [files(ii).folder,'/', files(ii).name];
    % open the file
    rsk = RSKopen(fin);
    %
    rsk = RSKreaddata(rsk);
    %  mat = RSK2MAT(rsk);
    serialNum = rsk.instruments.serialID;
    time      = rsk.data.tstamp;
    pres      = rsk.data.values;    
    % time in water
% $$$     figure, plot(time,pres)
% $$$     if ~useStartStop | ~exist('startStop','var')
% $$$         disp(['select a start and end time when instrument was in the water'])
% $$$         startStop = ginput(2);
% $$$     end
% $$$     it        = find(time>=startStop(1,1) &...
% $$$                      time<=startStop(2,1));
% $$$     pres      = pres(it);
% $$$     time      = time(it);
% $$$     % time in air
% $$$     disp(['select a start and end time when instrument was in the air'])
% $$$     startStop = ginput(2);
% $$$     it0       = find(time>=startStop(1,1) &...
% $$$                      time<=startStop(2,1));
% $$$     time0     = mean(rsk.data.tstamp(it0));
% $$$     pres0     = mean(rsk.data.values(it0));
    time0     = mean(rsk.data.tstamp);
    pres0     = mean(rsk.data.values);
    %
    % save data to structure
    data(ii).serialNum = serialNum;
    data(ii).time      = time;
    data(ii).pres      = pres;
    data(ii).time0     = time0;
    data(ii).pres0     = pres0;
end

save('../rbr_cal_20240207.mat','data')

