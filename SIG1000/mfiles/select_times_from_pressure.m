 % load the first file, find the deploy times and atmospheric pressure periods
load([rootDIR,filesep,files(1).name],'Data','Config','Descriptions')
%
Pres = Data.Burst_Pressure;
Time = Data.Burst_Time;
if isempty(atmosphTime)
    close all,
    figure, plot(Pres), title('select times when instrument was in air'),drawnow
    inAir       = ginput(2);
    atmosphTime = Time(round(inAir(:,1)))';
end
%
% log the pre-deployment atmosperic pressure time
inAir = (Time>=atmosphTime(1,1) & Time<=atmosphTime(1,2));
ATM_Time     = mean(Time(inAir));
ATM_Pressure = mean(Pres(inAir));
%
%
if isempty(deployTime)
    close all,
    figure, plot(Pres), title('select earliest time when the instrument was in water'),drawnow
    inWater     = ginput(1);
    TimeInWater = Time(round(inWater(:,1)));
    deployTime  = (TimeInWater*24-mod(TimeInWater*24,1)+1)/24;
    fprintf(' deployTime = %s \n',datestr(deployTime))
end
%
%
load([rootDIR,filesep,fRoot,num2str(Nf),'.mat'],'Data','Config')
Pres = Data.Burst_Pressure;
Time = Data.Burst_Time;
iter    = 0;
if isempty(recoverTime) | isempty(atmosphTime)
    inWater = [];
    while isempty(inWater)
        load([rootDIR,filesep,fRoot,num2str(Nf-iter),'.mat'],'Data','Config')
        Pres = Data.Burst_Pressure;
        Time = Data.Burst_Time;
        %
        close all,
        figure, plot(Pres), title('select latest time when the instrument was in water, if none press enter'),drawnow
        inWater      = ginput(1);
        if ~isempty(inWater) & isempty(recoverTime)
            TimeInWater  = Time(round(inWater(:,1)));
            recoverTime  = (TimeInWater*24-mod(TimeInWater*24,1))/24;
            fprintf(' recoverTime = %s \n',datestr(recoverTime))    
        end
        %
        close all,
        figure, plot(Pres), title('select times when instrument was in air, if none press enter'),drawnow
        inAir     = ginput(2);
        if ~isempty(inAir)
            atmosphTime = cat(1,atmosphTime,Time(round(inAir(:,1)))');
        end
        iter = iter+1;
    end
end

inAir   = (Time>=atmosphTime(end,1) & Time<=atmosphTime(end,2));    
while sum(inAir)==0 & ~isempty([inAir])
    iter = iter+1;
    load([rootDIR,filesep,fRoot,num2str(Nf-iter),'.mat'],'Data','Config')
    Pres = Data.Burst_Pressure;
    Time = Data.Burst_Time;
    inAir   = (Time>=atmosphTime(end,1) & Time<=atmosphTime(end,2));
% $$$     if sum(inAir)~=0 & ~isempty([inAir])
% $$$         atmosphTime = cat(1,atmosphTime,Time(round(inAir(:,1)))');
% $$$     end
end


close all

ATM_Time        = cat(1,ATM_Time    , mean(Time(inAir)));
ATM_Pressure    = cat(1,ATM_Pressure, mean(Pres(inAir)));

fprintf(' atmosphTime = \n')
datestr(atmosphTime)
