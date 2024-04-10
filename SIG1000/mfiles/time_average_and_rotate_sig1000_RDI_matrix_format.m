%
load([outDir,filesep,filePrefix,'config.mat'])
fs = double(Config.Burst_SamplingRate);
Nc = Config.Burst_NCells;
%
files    = dir([outDir,filesep,filePrefix,'*.mat']);
fNameCell=extractfield(files,'name');
files    = files(~contains(fNameCell,'config') & ~contains(fNameCell,'min.mat'));
Nf       = length(files);
%
fprintf('\n \n')
%
% define filter parameters
fc = 1/dtAvg;% frequency cutoff
Ns = fs/fc;% filter half-width
if isodd(Ns), Ns=Ns+1; end
Fw = 2*Ns+1;% filter width
F  = hanning(Ns+1);% Hann-function filter
F  = F./sum(F);% normalize
%
for ii= 1:Nf
    %
    % load the raw data
    inFileName = sprintf([filePrefix,'%03d.mat'],ii);
    fin  = [outDir,inFileName];
    fprintf(['loading file:   %s \n'],inFileName)
    in = load(fin,'Velocity_East','Velocity_North','Velocity_Up','Velocity_Error','Correlation_Minimum','Amplitude_Minimum','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','Temperature','bin_mab','bin_mab_Echo','Echo1','Echo2');
    %
    %
    if ii==1
        % no padding
        tP   = [];
        qcP  = [];
        vb1p = [];
        vb2p = [];
        vb3p = [];
        vb4p = [];
        e1p  = [];
        e2p  = [];
        hP   = [];
        pP   = [];
        rP   = [];
        Pp   = [];
        Tp   = [];
    end
    %
    np = length(tP);
    nt = length(in.Time);
    % fractional number of full filter segements
    ns = (nt+np)/Ns;
    % number of whole segments
    N  = floor( ns );
    % remainder
    ns = rem(ns,1);
    %
    % average Iburst and burst Temp/Pres
    in.Pressure    = nanmean(in.Pressure,1);
    in.Temperature = nanmean(in.Temperature,1);
    %
    %
    % perform time-convolution (filling nan's in the mask), conv2 x 3 (nans=0, qcFlag, ones)
    on  = [qcP, in.qcFlag];
    o   = ones(1,size(on,2));
    %
    % convolve NAN matrix for normalization
    norm0   = conv2(o,F,'same');
    avgFlag = conv2(on,F,'same')./norm0;
    % throw out regions with <%25 coverage
    norm = avgFlag;
    norm(norm<=0.25)=inf;
    %
    % convolve signals
    t1    = [tP, in.Time];
    vb1   = conv2(1, F, [vb1p, in.Velocity_East ].*on,'same')./norm;
    vb2   = conv2(1, F, [vb2p, in.Velocity_North].*on,'same')./norm;
    vb3   = conv2(1, F, [vb3p, in.Velocity_Up   ].*on,'same')./norm;
    vb4   = conv2(1, F, [vb4p, in.Velocity_Error].*on,'same')./norm;
    echo1 = conv2(1, F, [e1p , in.Echo1]       ,'same')./norm0;
    echo2 = conv2(1, F, [e2p , in.Echo2]       ,'same')./norm0;    
    head  = conv2(1, F, [hP  , in.Heading ]    ,'same')./norm0;
    pitch = conv2(1, F, [pP  , in.Pitch   ]    ,'same')./norm0;
    roll  = conv2(1, F, [rP  , in.Roll    ]    ,'same')./norm0;
    P     = conv2(1, F, [Pp  , in.Pressure]    ,'same')./norm0;
    T     = conv2(1, F, [Tp  , in.Temperature] ,'same')./norm0;    
    %
    %
    avg = struct('Time',t1(Ns:Ns:Ns*(N-1)),'Velocity_East',vb1(:,Ns:Ns:Ns*(N-1)),'Velocity_North',vb2(:,Ns:Ns:Ns*(N-1)),'Velocity_Up',vb3(:,Ns:Ns:Ns*(N-1)),'Velocity_Error',vb4(:,Ns:Ns:Ns*(N-1)),'Heading',head(:,Ns:Ns:Ns*(N-1)),'Pitch',pitch(:,Ns:Ns:Ns*(N-1)),'Roll',pitch(:,Ns:Ns:Ns*(N-1)),'Pressure',P(Ns:Ns:Ns*(N-1)),'Temperature',T(Ns:Ns:Ns*(N-1)),'qcFlag',avgFlag(:,Ns:Ns:Ns*(N-1)),'HeadingOffset',in.HeadingOffset,'bin_mab',in.bin_mab,'bin_mab_Echo',in.bin_mab_Echo,'Echo1',echo1(:,Ns:Ns:Ns*(N-1)),'Echo2',echo2(:,Ns:Ns:Ns*(N-1)));    
    %
    % 
    % need to rotate from magnetic to true?
% $$$     [avg,Config2,beam2xyz] = beam2earth_sig1000_DG(avg,Config,'',in.HeadingOffset);
% $$$     avg.beam2xyz=beam2xyz;
    %
    if ii==1
        out = avg;
    else
        vars = fields(out);
        for jj = 1:length(vars);
            var = vars{jj};
            if ismember(var,{'HeadingOffset','bin_mab','bin_mab_Echo','beam2xyz'})
                continue
            end
            eval(['out.',var,'= cat(2,out.',var,',avg.',var,');'])
        end
    end
    %
    % get pad data for next file
    tP   = in.Time       (:,nt-Ns*(1+ns)+1:nt);
    qcP  = in.qcFlag     (:,nt-Ns*(1+ns)+1:nt);
    vb1p = in.Velocity_East    (:,nt-Ns*(1+ns)+1:nt);
    vb2p = in.Velocity_North   (:,nt-Ns*(1+ns)+1:nt);
    vb3p = in.Velocity_Up      (:,nt-Ns*(1+ns)+1:nt);
    vb4p = in.Velocity_Error   (:,nt-Ns*(1+ns)+1:nt);
    e1p  = in.Echo1      (:,nt-Ns*(1+ns)+1:nt);
    e2p  = in.Echo2      (:,nt-Ns*(1+ns)+1:nt);    
    hP   = in.Heading    (nt-Ns*(1+ns)+1:nt);
    pP   = in.Pitch      (nt-Ns*(1+ns)+1:nt);
    rP   = in.Roll       (nt-Ns*(1+ns)+1:nt);
    Tp   = in.Temperature(nt-Ns*(1+ns)+1:nt);
    Pp   = in.Pressure   (nt-Ns*(1+ns)+1:nt);
    %
end
fprintf(['saving: %s'],[L0dir,L0FRoot])
if ~exist(L0dir,'dir')
    eval(['!mkdir -p ',L0dir])
end
save([L0dir,L0FRoot],'-struct','out')
disp('done!')