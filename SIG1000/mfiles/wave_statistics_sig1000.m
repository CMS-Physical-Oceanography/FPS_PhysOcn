%
dtBurst = 1800;% seconds
dtEns   = 512 ;% seconds
%
%
% get sampling information from config file
load([outRoot,L0FRoot,'config.mat'])
fs = double(Config.Burst_SamplingRate);
Nc = Config.Burst_NCells;
%
% sample time step
dt = 1/fs;
%
% number of samples to average
Na = dtBurst*fs;
% number of samples in each ensemble
Ne    = dtEns*fs;
olap  = 2/3;
chnks = (Na-Ne*olap-1)/(Ne*(1-olap));
%
% get structure with all files in archive
files = dir([outRoot,L0FRoot,'*.mat']);
fNameCell=extractfield(files,'name');
files = files(~contains(fNameCell,'config') & ~contains(fNameCell,'min.mat'));
Nf    = length(files);
%
% initialize counter...
fprintf('\n \n')
ensNum = 1;
out = struct();
for ii= 1:Nf
    %
    % load the raw data
    inFileName = sprintf([L0FRoot,'%03d.mat'],ii);
    fin  = [outRoot,inFileName];
    fprintf(['loading file:   %s \n'],inFileName)
% $$$     in = load(fin,'VelEast','VelNorth','VelUp1','VelUp2','VelBeam5','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','mab');
    in = load(fin,'VelEast','VelNorth','VelUp','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','mab');
    %
    if ii==1
        mab = in.mab;
    end
    %
    % get number of ensembles in current file:
    nt = length(in.Time);
    N = floor(nt/Na);
    % 
    % loop over ensembles
    for jj = 1:N
        %
        % check to see which bins are in the water
        QC   = in.qcFlag((jj-1)*Na+1:Na*jj,:);
        fGood= sum(QC,1)/Na;
        % only use bins in water >99 % of the time
        bins = find(fGood>0.99);
        Nb   = length(bins);
        %
        % mean time to nearest minute
        t = in.Time((jj-1)*Na+1:Na*jj);
        tavg = round(mean(t)*1440)/1440;
        %
        out.Time_wave(ensNum,1) = tavg;
        if Nb<2
            disp('not enough good data')
            out.Hs(ensNum,1) = nan;
            out.Tm(ensNum,1) = nan;
            out.Suu(ensNum,:)= nan;
            out.Svv(ensNum,:)= nan;
            out.Spp(ensNum,:)= nan;
            out.Spu(ensNum,:)= nan;
            out.Spp(ensNum,:)= nan;
            out.Spv(ensNum,:)= nan;
% $$$             out.cP2eta(ensNum,:)  = nan;        
% $$$             out.cU2eta(ensNum,:,:)  = nan;
% $$$             out.cU2srf(ensNum,:,:)  = nan;
            out.mdir(ensNum,1:2)  = nan;
            out.mspread(ensNum,1:2)=nan;
            out.mSxx(ensNum) = nan;
            out.mSxy(ensNum) = nan;
            out.mSyy(ensNum) = nan;
            out.Z2(ensNum,:) = nan;
            ensNum = ensNum+1;
            continue
        end
        %
        % allocate current variables: (note, bad QC are zeros not nans, and not interpolated)
        U = in.VelEast ((jj-1)*Na+1:Na*jj,bins);
        V = in.VelNorth((jj-1)*Na+1:Na*jj,bins);
        W = in.VelUp   ((jj-1)*Na+1:Na*jj,bins);
        P = in.Pressure((jj-1)*Na+1:Na*jj,2);
        %
        % depth of pressure sensor
        dpthP= mean(P,1);% not accounting for sensor height above bed
        % assuming sensor is at bed level        
        dpth = dpthP;
        % depth of bins
        dpthU= dpthP-mab(bins);
        %
        % estimate bulk wave statistics on 10 min chuncks, no overlap
        [Suu,fq]=welch_method(U,dt,chnks,olap);
        [Svv,fq]=welch_method(V,dt,chnks,olap);
        [Spp,fq]=welch_method(P,dt,chnks,olap);
        [Suv,fq]=welch_cospec(U,V,dt,chnks,olap);        
        [Spu,fq]=welch_cospec(repmat(P,1,Nb),U,dt,chnks,olap);
        [Spv,fq]=welch_cospec(repmat(P,1,Nb),V,dt,chnks,olap);
        %
        % note: still has zero frequency!
        fq = fq(2:end);
        Suu= Suu(2:end,:);
        Svv= Svv(2:end,:);
        Spp= Spp(2:end,:);
        Suv= Suv(2:end,:);
        Spu= Spu(2:end,:);
        Spv= Spv(2:end,:);
        %
        % depth correction
        om = 2*pi*fq;
        lom=length(om);
        %
        k=wavenumber(om.',dpth); % these are the radian wave numbers (rad/m)
        g=9.81;
        %
        % convert pressure to surface elevation
        cP     =  cosh(k.*dpth)./ cosh(k.*(dpth-dpthP)); % SePP=Spp.*cP.^2
        % convert velocity to surface elevation
        cU2eta = (om./(g*k)).*cosh(k.*dpth)./cosh(k.*(dpth-dpthU));% SeUU = SeUU.*cU.^2        
        % convert velocity to surface velocity
        cU2srf =              cosh(k.*dpth)./cosh(k.*(dpth-dpthU));% SeUU = SeUU.*cU.^2        
        %
        % estimate max resolvable frequency where depth bin is below wave-length/(pi)
        fmax = min(fq(Ne/2-sum(dpthU>=1./(k)))');
        % nan surface transformed fields above fmax
        Suu(fq>fmax,:)=nan;
        Svv(fq>fmax,:)=nan;
        Suv(fq>fmax,:)=nan;
        Spp(fq>fmax(1))=nan;
        Spu(fq>fmax,:)=nan;
        Spv(fq>fmax,:)=nan;
        % surface velocity spectra
        SUU = Suu.*cU2srf.^2;
        SVV = Svv.*cU2srf.^2;
        SUV = Suv.*cU2srf.^2;
        SPU = repmat(cP,1,Nb).*Spu.*cU2srf;
        SPV = repmat(cP,1,Nb).*Spv.*cU2srf;        
        %
        % map to surface elevation spectra
        SeUU = Suu.*cU2eta.^2;
        SeVV = Svv.*cU2eta.^2;
        SeUV = Suv.*cU2eta.^2;
        SePP = Spp.*cP.^2;
        SePU = repmat(cP,1,Nb).*Spu.*cU2eta;
        SePV = repmat(cP,1,Nb).*Spv.*cU2eta;
        %
        %
        if hour(tavg)==12
            fig1 = figure;
            plt = semilogy(fq,mean(SeUU+SeVV,2),'-r',fq,SePP,'--k','linewidth',3);
            set(gca,'xlim',[1/(30*60) 1/3])
            grid on
            xlabel('cycles per second','interpreter','latex')
            ylabel('m$^2$/Hz','interpreter','latex')
            title(datestr(tavg),'interpreter','latex')
            h = legend(plt,'$S_{uu}+S_{vv}$','$S_{\eta\eta}$');
            set(h,'location','northeast','interpreter','latex')
            figname = sprintf('../figures/%s_spectra_%s.pdf',L1FRoot,datestr(tavg,'mmddHHMM'));
            exportgraphics(fig1,figname)
            close(fig1)
        end
        %
        % surface corrected total velocity spectra
        SUtot = SeUU+SeVV;
        Z2 = nanmean(SePP./SUtot,2);
        %
        %
% $$$         % create fake (P,U,V), rotated by -70deg=160 nautical, or +70deg=20nautical
% $$$         Nb = 1;
% $$$         a  = 1;% 1 meter
% $$$         ang= 70*pi/180+pi;% direction from--> direction to
% $$$         P  = a*cos(2*pi/10.*(t-t(1))*86400);
% $$$         U  = P.*cos(ang) ;
% $$$         V  = P.*sin(ang);
% $$$         %
% $$$         % estimate bulk wave statistics on 10 min chuncks, no overlap
% $$$         [Suu,fq]=welch_method(U,dt,chnks,olap);
% $$$         [Svv,fq]=welch_method(V,dt,chnks,olap);
% $$$         [Spp,fq]=welch_method(P,dt,chnks,olap);
% $$$         [Suv,fq]=welch_cospec(U,V,dt,chnks,olap);        
% $$$         [Spu,fq]=welch_cospec(repmat(P,1,Nb),U,dt,chnks,olap);
% $$$         [Spv,fq]=welch_cospec(repmat(P,1,Nb),V,dt,chnks,olap);
% $$$         %
% $$$         % note: still has zero frequency!
% $$$         fq  = fq(2:end);
% $$$         SUU = Suu(2:end,:);
% $$$         SVV = Svv(2:end,:);
% $$$         SePP= Spp(2:end,:);
% $$$         SUV = Suv(2:end,:);
% $$$         SPU = Spu(2:end,:);
% $$$         SPV = Spv(2:end,:);
% $$$         %
        % estimate bulk stats
        % I think there is a sign issue here
        coPU =  SPU;
        coPV =  SPV;
        coUV =  SUV;
        r2d = 180/pi;
        %
        % get a1
        a1      = coPU ./ sqrt( SePP .* ( SUU + SVV ) );
        b1      = coPV ./ sqrt( SePP .* ( SUU + SVV ) );
        dir1    = r2d* ( atan2(b1,a1) );          
        spread1 = r2d* ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );
        %
        % average over wind-wave band
        df= fq(2)-fq(1);
        I = find(fq>=1/20 & fq<=min(1/4,min(fmax)));
        m0 = nansum(SePP(I)*df);
        m1 = nansum(fq(I).*SePP(I)*df);        
        ma1= nansum(a1(I,:).*SePP(I)*df,1)/m0;
        mb1= nansum(b1(I,:).*SePP(I)*df,1)/m0;
        mdir1=rem(90+180-r2d*atan2(mb1,ma1),360);% nautical 
        mspread1 = r2d*sqrt(abs(2*(1-(ma1.*cos(mdir1/r2d) + mb1.*sin(mdir1/r2d)))));
        %
        % get a2 b2
        a2 = (SUU - SVV) ./ (SUU + SVV);
        b2 = 2 .* coUV   ./ (SUU + SVV);
        spread2 = r2d*sqrt(abs(0.5-0.5*(a2.*cos(2.*dir1/r2d)+b2.*sin(2.*dir1/r2d))));
        %
        ma2      = nansum(a2(I,:).*SePP(I)*df,1)/m0;
        mb2      = nansum(b2(I,:).*SePP(I)*df,1)/m0;
        dir2     = r2d/2*atan2(b2,a2);
        mdir2    = 90-(r2d/2*atan2(mb2,ma2));% nautical
        mspread2 = r2d*sqrt(abs(0.5-0.5*(ma2.*cos(2.*mdir1/r2d)+mb2.*sin(2.*mdir1/r2d))));
        %
        % radiation stresses
        % Radiation Stress Estimates
        C = om./k;
        Cg = 0.5*(g*tanh(k.*dpth)+g*k.*dpth.*(sech(k.*dpth).^2))./sqrt(g*k.*tanh(k.*dpth));
        %
        % cartesian
        Sxx = ( (1.5 + 0.5*a2) .* (Cg./C) - 0.5 ) .* SePP;
        Syy = ( (1.5 - 0.5*a2) .* (Cg./C) - 0.5 ) .* SePP;
        Sxy = 0.5*b2 .* (Cg./C) .* SePP;
        %
        mSxx = sum( SePP(I).*Sxx(I,:).*df )./m0;
        mSyy = sum( SePP(I).*Syy(I,:).*df )./m0;
        mSxy = sum( SePP(I).*Sxy(I,:).*df )./m0;                
        %
        %
        Hs = 4*sqrt(nansum(df*SePP(I)));
        Tm = m0/m1;
        %
        % log the hourly averaged
        if ii==1 & jj==1
            out.frequency  = fq';
            out.wavenumber = k';            
        end
        out.Hs     (ensNum,1)   = Hs;
        out.Tm     (ensNum,1)   = Tm;
        out.depth  (ensNum,1)   = dpth;
        out.Spp    (ensNum,:)   = SePP;
        out.Suu    (ensNum,:)   = nanmean(SUU,2);
        out.Svv    (ensNum,:)   = nanmean(SVV,2);
        out.Spu    (ensNum,:)   = nanmean(SPU,2);
        out.Spv    (ensNum,:)   = nanmean(SPV,2);
        out.mdir   (ensNum,1:2) = [nanmean(mdir1)   ,nanmean(mdir2)];
        out.mspread(ensNum,1:2) = [nanmean(mspread1),nanmean(mspread2)];
        out.mSxx   (ensNum)     = nanmean(mSxx);
        out.mSxy   (ensNum)     = nanmean(mSxy);
        out.mSyy   (ensNum)     = nanmean(mSyy);
        out.Z2    (ensNum,:)    = Z2;
        ensNum = ensNum+1;
    end
    %
end
fprintf(['saving: %s'],[outRoot,L1FRoot])
save([outRoot,L1FRoot,'.mat'],'-append','-struct','out')
disp('done!')