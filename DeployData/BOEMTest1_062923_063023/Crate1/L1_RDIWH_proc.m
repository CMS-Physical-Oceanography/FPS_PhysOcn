%% Section 3 Finding the time the ADCP is in and out the water
% use depth field first to try and determine pre and post deployment data

for i = 1:length(ADCP.cur)
       clear bad_depth
       if length(unique(ADCP.cur(i).depth))>2 && mean(ADCP.cur(i).depth) > 12 %Checking to make sure pressure sensor is recording and in the range of what we expect
           
           bad_depth = find(ADCP.cur(i).depth < 9.5); %used to index values that are assumed to be above the surface. All values where depth is less than 5 are deleted if they meet the if statement above

           ADCP.cur(i).mtime(bad_depth)=[];
           ADCP.cur(i).pitch(bad_depth)=[];
           ADCP.cur(i).roll(bad_depth)=[];
           ADCP.cur(i).heading(bad_depth)=[];
           ADCP.cur(i).depth(bad_depth)=[];
           ADCP.cur(i).temperature(bad_depth)=[];
           ADCP.cur(i).salinity(bad_depth)=[];
           ADCP.cur(i).pressure(bad_depth)=[];
           ADCP.cur(i).east_vel(:,bad_depth)=[];
           ADCP.cur(i).north_vel(:,bad_depth)=[];
           ADCP.cur(i).vert_vel(:,bad_depth)=[];
           ADCP.cur(i).error_vel(:,bad_depth)=[];
           ADCP.cur(i).corr(:,:,bad_depth)=[];
           ADCP.cur(i).intens(:,:,bad_depth)=[];
           ADCP.cur(i).perc_good(:,:,bad_depth)=[];
       else
       end
    clear bad_depth
end

%The echo amplitude is now used to determine pre and post deployment data.
%The echo amplitude (intensity) is expected to be lower in air than in
%water. With a 0.5 meter bin size, Bin 25 is near the surface and has a high intensity, if it is low we
%assume that the instrument is out of the water.
for i = 1:length(ADCP.cur)
       clear bad_echo
       
           bad_echo = find(ADCP.cur(i).intens(25,1,[1:15 end-15:end]) < 90); 
           %intens(#, , ) needs to be set at the surface intensities where it begin to spike from the surface
           ADCP.cur(i).mtime(bad_echo)=[];
           ADCP.cur(i).pitch(bad_echo)=[];
           ADCP.cur(i).roll(bad_echo)=[];
           ADCP.cur(i).heading(bad_echo)=[];
           ADCP.cur(i).depth(bad_echo)=[];
           ADCP.cur(i).temperature(bad_echo)=[];
           ADCP.cur(i).salinity(bad_echo)=[];
           ADCP.cur(i).pressure(bad_echo)=[];
           ADCP.cur(i).east_vel(:,bad_echo)=[];
           ADCP.cur(i).north_vel(:,bad_echo)=[];
           ADCP.cur(i).vert_vel(:,bad_echo)=[];
           ADCP.cur(i).error_vel(:,bad_echo)=[];
           ADCP.cur(i).corr(:,:,bad_echo)=[];
           ADCP.cur(i).intens(:,:,bad_echo)=[];
           ADCP.cur(i).perc_good(:,:,bad_echo)=[];
end
disp('Section 3 complete')
%% Section 4 Manual input to determine pre and post deployment
%the instruments heading is plotted and the user will determine the moments
%for pre and post data

figure
for i = 1:length(ADCP.cur)
    clear post_data pre_data samp_int xx yy
    
    subplot(2,1,1) %top plot is of the first 50 ensembles heading that remain after testing the depth and echo amplitude
    plot(ADCP.cur(i).mtime,ADCP.cur(i).heading,'Linewidth',2)
    xlim([ADCP.cur(i).mtime(1) ADCP.cur(i).mtime(50)])
    ylabel('degrees from north CW')
    title('Heading for First 50 Ensembles')
    xlabel('ensemble number')
    grid on

    subplot(2,1,2)% bottom plot is of the last 50 ensembles of deployment
    plot(ADCP.cur(i).mtime,ADCP.cur(i).heading,'Linewidth',2)
    xlim([ADCP.cur(i).mtime(end-50) ADCP.cur(i).mtime(end)])
    ylabel('degrees from north CW')
    title('Heading for Last 50 Ensembles')
    xlabel('ensemble number from end')
    grid on
    
    [xx,yy] = ginput(2); %first click needs to be on the top graph AFTER the heading becomes stable. second click is on the bottom graph BEFORE the heading becomes unstable
    %only the x coordinate will matter
    
    samp_int = ADCP.cur(i).mtime(10)-ADCP.cur(i).mtime(9); %determines the sample interval for the instruments deployment
    xx=[xx(1)+2*samp_int xx(2)-2*samp_int]; %takes the time determined from the user's input on heading graph and adds two sampling intervals before the first click and after the second click
    %the two sampling intervals are added to ensure the bad ensemble is
    %captured and deleted
    
    pre_data = find(ADCP.cur(i).mtime < xx(1)); %pre deployment data 
    post_data = find(ADCP.cur(i).mtime > xx(2)); %post deployment data
    
    %post deployment data is deleted first to not shift the indexes one
    %way. If the pre deployment data was deleted first the indexs for the
    %post deployment data would all be off
     ADCP.cur(i).mtime(post_data)=[];
     ADCP.cur(i).pitch(post_data)=[];
     ADCP.cur(i).roll(post_data)=[];
     ADCP.cur(i).heading(post_data)=[];
     ADCP.cur(i).depth(post_data)=[];
     ADCP.cur(i).temperature(post_data)=[];
     ADCP.cur(i).salinity(post_data)=[];
     ADCP.cur(i).pressure(post_data)=[];
     ADCP.cur(i).east_vel(:,post_data)=[];
     ADCP.cur(i).north_vel(:,post_data)=[];
     ADCP.cur(i).vert_vel(:,post_data)=[];
     ADCP.cur(i).error_vel(:,post_data)=[];
     ADCP.cur(i).corr(:,:,post_data)=[];
     ADCP.cur(i).intens(:,:,post_data)=[];
     ADCP.cur(i).perc_good(:,:,post_data)=[];
    
     ADCP.cur(i).mtime(pre_data)=[];
     ADCP.cur(i).pitch(pre_data)=[];
     ADCP.cur(i).roll(pre_data)=[];
     ADCP.cur(i).heading(pre_data)=[];
     ADCP.cur(i).depth(pre_data)=[];
     ADCP.cur(i).temperature(pre_data)=[];
     ADCP.cur(i).salinity(pre_data)=[];
     ADCP.cur(i).pressure(pre_data)=[];
     ADCP.cur(i).east_vel(:,pre_data)=[];
     ADCP.cur(i).north_vel(:,pre_data)=[];
     ADCP.cur(i).vert_vel(:,pre_data)=[];
     ADCP.cur(i).error_vel(:,pre_data)=[];
     ADCP.cur(i).corr(:,:,pre_data)=[];
     ADCP.cur(i).intens(:,:,pre_data)=[];
     ADCP.cur(i).perc_good(:,:,pre_data)=[];
     
end 
disp('Section 4 complete')
%% Section 5 Get rid of bad pressure
%a graph of pressure will be displayed and the user will be asked to input
%a 1 if the pressure is good or a 0 if the pressure is bad
%
%If the pressure reading is wacky it is deleted 
figure
for i = 1:length(ADCP.cur)
    plot(ADCP.cur(i).mtime,ADCP.cur(i).pressure)
    datetick
    xlim([ADCP.cur(i).mtime(1) ADCP.cur(i).mtime(end)])
%     ylim([2.4e4 3e4])
%     yticks([2.4e4:1e4:3e4])
%     yticklabels({'24000','25000','26000','27000','29000','30000'})
    ylabel('pressure (dBars)')
    xlabel('date')
    title('Time Series of Presure')
    grid on
    
    answer = inputdlg('Is Pressure Data Good? (1 for yes 0 for no)');

    if str2double(answer{1}) == 0
       ADCP.cur(i).pressure = [];
    else
    end
end       
disp('Section 5 complete')
%% %% Section 7 Quality Check on Section 6
% if test image shows something funky run this section
dep_num = inputdlg('Which deployment do you wish to trim? Enter 1 if there is only 1 deployment');
i = str2double(dep_num{1}); % make i equal to the deployment you wish to edit. If you are dealing with only 1 deployment i = 1.

beg_bad = inputdlg('How many ensembles need to be trimmed from begining of time series?');
end_bad = inputdlg('How many ensembles need to be trimmed from end of time series?');

beg_bad = str2double(beg_bad{1});
end_bad = str2double(end_bad{1});

if beg_bad ~= 0
    trim_start = [1:beg_bad]; %manually input the number of ensembles that appear to have bad data and they will be deleted from record.
    
    %the following if-else loop is used because if there is pressure data the
    %if statement is used to deleted the selected ensembles. 
    %%%
    % However, if there is no pressure data the else statement deleted the bad
    % ensembles. This is needed not to have an error since bad pressure data
    % was already deleted 
    if isempty(ADCP.cur(i).pressure) == 0
    
     ADCP.cur(i).mtime(trim_start)=[];
     ADCP.cur(i).pitch(trim_start)=[];
     ADCP.cur(i).roll(trim_start)=[];
     ADCP.cur(i).heading(trim_start)=[];
     ADCP.cur(i).depth(trim_start)=[];
     ADCP.cur(i).temperature(trim_start)=[];
     ADCP.cur(i).salinity(trim_start)=[];
     ADCP.cur(i).pressure(trim_start)=[];
     ADCP.cur(i).east_vel(:,trim_start)=[];
     ADCP.cur(i).north_vel(:,trim_start)=[];
     ADCP.cur(i).vert_vel(:,trim_start)=[];
     ADCP.cur(i).error_vel(:,trim_start)=[];
     ADCP.cur(i).corr(:,:,trim_start)=[];
     ADCP.cur(i).intens(:,:,trim_start)=[];
     ADCP.cur(i).perc_good(:,:,trim_start)=[];
    else
     ADCP.cur(i).mtime(trim_start)=[];
     ADCP.cur(i).pitch(trim_start)=[];
     ADCP.cur(i).roll(trim_start)=[];
     ADCP.cur(i).heading(trim_start)=[];
     ADCP.cur(i).depth(trim_start)=[];
     ADCP.cur(i).temperature(trim_start)=[];
     ADCP.cur(i).salinity(trim_start)=[];
     ADCP.cur(i).east_vel(:,trim_start)=[];
     ADCP.cur(i).north_vel(:,trim_start)=[];
     ADCP.cur(i).vert_vel(:,trim_start)=[];
     ADCP.cur(i).error_vel(:,trim_start)=[];
     ADCP.cur(i).corr(:,:,trim_start)=[];
     ADCP.cur(i).intens(:,:,trim_start)=[];
     ADCP.cur(i).perc_good(:,:,trim_start)=[];
    end
else
end

if end_bad ~= 0
    trim_end = [length(ADCP.cur.mtime)-end_bad:length(ADCP.cur.mtime)];
    
    %the following if-else loop is used because if there is pressure data the
    %if statement is used to deleted the selected ensembles. 
    %%%
    % However, if there is no pressure data the else statement deleted the bad
    % ensembles. This is needed not to have an error since bad pressure data
    % was already deleted 
    if isempty(ADCP.cur(i).pressure) == 0
    
     ADCP.cur(i).mtime(trim_end)=[];
     ADCP.cur(i).pitch(trim_end)=[];
     ADCP.cur(i).roll(trim_end)=[];
     ADCP.cur(i).heading(trim_end)=[];
     ADCP.cur(i).depth(trim_end)=[];
     ADCP.cur(i).temperature(trim_end)=[];
     ADCP.cur(i).salinity(trim_end)=[];
     ADCP.cur(i).pressure(trim_end)=[];
     ADCP.cur(i).east_vel(:,trim_end)=[];
     ADCP.cur(i).north_vel(:,trim_end)=[];
     ADCP.cur(i).vert_vel(:,trim_end)=[];
     ADCP.cur(i).error_vel(:,trim_end)=[];
     ADCP.cur(i).corr(:,:,trim_end)=[];
     ADCP.cur(i).intens(:,:,trim_end)=[];
     ADCP.cur(i).perc_good(:,:,trim_end)=[];
    else
     ADCP.cur(i).mtime(trim_end)=[];
     ADCP.cur(i).pitch(trim_end)=[];
     ADCP.cur(i).roll(trim_end)=[];
     ADCP.cur(i).heading(trim_end)=[];
     ADCP.cur(i).depth(trim_end)=[];
     ADCP.cur(i).temperature(trim_end)=[];
     ADCP.cur(i).salinity(trim_end)=[];
     ADCP.cur(i).east_vel(:,trim_end)=[];
     ADCP.cur(i).north_vel(:,trim_end)=[];
     ADCP.cur(i).vert_vel(:,trim_end)=[];
     ADCP.cur(i).error_vel(:,trim_end)=[];
     ADCP.cur(i).corr(:,:,trim_end)=[];
     ADCP.cur(i).intens(:,:,trim_end)=[];
     ADCP.cur(i).perc_good(:,:,trim_end)=[];
    end
else 
end

%the following if-else loop is used because if there is pressure data the
%if statement is used to deleted the selected ensembles. 
%%%
% However, if there is no pressure data the else statement deleted the bad
% ensembles. This is needed not to have an error since bad pressure data
% was already deleted
disp('Section 7 complete')          
%% Section 8 Applying QC procedures
% north values
% velocity error test with a threshold of 0.05 cm/s
for zz=1:length(ADCP.cur)
    clear col dif ind row xx yy 
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
disp('Section 8 complete')
%% Section 9 East values
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
disp('Section 9 complete')
%% %% Section 11 Calculating 4 beam average Echo Amplitude
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
for zz = 1:length(ADCP.cur)
    % create a set of row indices
    inds = repmat([1:ADCP.cur(zz).config.n_cells]',1,length(ADCP.cur(zz).mtime));
    % get first difference in rows
    dum0 = ADCP.cur(zz).EA_avg(1:ADCP.cur(zz).config.n_cells-1,:);
    dum1 = ADCP.cur(zz).EA_avg(2:ADCP.cur(zz).config.n_cells  ,:);
    dif  = dum1-dum0;
    % set the indices of any row with abs(dif)<20 to infinity
    inds(abs(dif)<20)=inf;
    c    = min(inds,[],1);
%     for xx=1:length(ADCP.cur(zz).mtime)
%         for yy=15:ADCP.cur(zz).config.n_cells-1 %only bins above the 15th bin
%             dif(xx,yy) = ADCP.cur(zz).EA_avg(yy,xx)-ADCP.cur(zz).EA_avg(yy+1,xx);
%            
%             if isempty(min(find(abs(dif(xx,:))>10))) == 0
%                 c(xx) = min(find(abs(dif(xx,:))>10));
%             else
%                 c(xx) = 25;
%             end
%         end
%     end
%  
%the config ranges is used to go from bin number to distance from ADCP
%(depth). then a mean across the entire deployment is taken
   ADCP.cur(zz).depth_ea = mean(abs(ADCP.cur(zz).config.ranges(c))); 
end
end

disp('Section 11 complete')
%% 
%% Section 13 Removing values that are above Sea surface 
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

    else
        surface_ea = find(abs(ADCP.cur.config.ranges) > ADCP.cur.depth_ea);
        surface_cut_ea = min(surface_ea);

        ADCP.cur.east_vel(surface_cut_ea:end,:) = nan;
        ADCP.cur.north_vel(surface_cut_ea:end,:) = nan;
    end
end
disp('Section 13 complete')




%% 
%======================================================================
% Selecting start time and stop time
datetime_vec=datetime(RBRtri.time, 'Format','yyyy-MM-dd HH:mm:ss','convertFrom','datenum');
datetime_strNend=datetime_vec([1:25 end-25:end]);
disp(datetime_strNend)

startTime=input('Enter Start Time yyyy-MM-dd HH:mm:ss: ','s')
endTime=input('Enter End Time yyyy-MM-dd HH:mm:ss: ','s')
start_ind=find(startTime==datetime_vec);
end_ind=find(datetime_vec==endTime);

%Modifying the lengths of Chlorophyll-a, FDOM, & Turbdidity to match time
RBRtri.time=datetime_vec(start_ind:end_ind);
RBRtri.chlorophyll_a.data=RBRtri.chlorophyll_a.data(start_ind:end_ind);
RBRtri.FDOM.data=RBRtri.FDOM.data(start_ind:end_ind);
RBRtri.turbidity.data=RBRtri.turbidity.data(start_ind:end_ind);

%If no input is made into the startTime and endTime, the function will plot
%the raw data using RSKplotdata
if isempty(startTime)==1 && isempty(endTime)==1
 RSKplotdata(RSK)

elseif  isempty(startTime)==0 && isempty(endTime)==0
  figure(); clf;
     subplot(311);
        plot(RBRtri.time,RBRtri.chlorophyll_a.data, 'LineWidth',1.5,'Color','blue');
        title('RBR Tridente Raw Data'); ylabel('phyll-a (ug/l)'); grid on;
        set(gca,'XTickLabel',[]);
    subplot(312);
        plot(RBRtri.time, RBRtri.FDOM.data, 'LineWidth',1.5,'Color','green');
        ylabel('FDOM (ppb)'); grid on; set(gca,'XTickLabel',[]);
    subplot(313);
        plot(RBRtri.time, RBRtri.turbidity.data, 'LineWidth',1.5,'Color','red');
        ylabel('turb (FTU)'); grid on;
end   



%% =======================================================================
% Selecting start time and stop time
datetime_vec=datetime(RBRtri.time, 'Format','yyyy-MM-dd HH:mm:ss','convertFrom','datenum');
datetime_strNend=datetime_vec([1:25 end-25:end]);
disp(datetime_strNend)

startTime=input('Enter Start Time yyyy-MM-dd HH:mm:ss: ','s')
endTime=input('Enter End Time yyyy-MM-dd HH:mm:ss: ','s')
start_ind=find(startTime==datetime_vec);
end_ind=find(datetime_vec==endTime);

%Modifying the structure to fit the desired start and end time
SB37.temperature.data=SB37.temperature.data(start_ind:end_ind);
SB37.conductivity.data=SB37.conductivity.data(start_ind:end_ind);
SB37.pressure.data=SB37.pressure.data(start_ind:end_ind);
SB37.time=SB37.time(start_ind:end_ind);

 if isempty(startTime)==1 && isempty(endTime)==1
  figure(); clf
    subplot(311)
    plot(SB37.time, SB37.temperature.data,'LineWidth',1.5,'Color','blue'); ylabel('Deg C');
    title('Raw SB37 Data'); grid on;
    subplot(312)
    plot(SB37.time, SB37.conductivity.data, 'LineWidth', 1.5, 'Color', 'red'); ylabel('S/m');
    grid on
    subplot(313)
    plot(SB37.time, SB37.pressure.data, 'LineWidth',1.5, 'Color', 'green'); ylabel('dB'); grid on
       
 else if isempty(starTime)==0 && isempty(endTime)==0
   figure(); clf;
    subplot(311)
    plot(time, temp,'LineWidth',1.5,'Color','blue'); ylabel('Deg C');
    title('Raw SB37 Data'); grid on;
    subplot(312)
    plot(time, conductivity, 'LineWidth', 1.5, 'Color', 'red'); ylabel('S/m');
    grid on
    subplot(313)
    plot(time, pres, 'LineWidth',1.5, 'Color', 'green'); ylabel('dB'); grid on