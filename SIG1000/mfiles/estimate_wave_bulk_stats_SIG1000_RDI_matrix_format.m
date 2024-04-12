 %
dtBurst = 1800;% seconds
dtEns   = 512 ;% seconds
rho0    = 1027.5;% kg/m^3
g       = 9.81;% m/s^2
%
%
% get sampling information from config file
load([outDir,filePrefix,'config.mat'])
fs = double(Config.Burst_SamplingRate);
Nc = Config.Burst_NCells;
%
% load L0-file with time-averages
currents = load([L0dir,filesep,L0FRoot,'.mat']);
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
files = dir([outDir,filePrefix,'*.mat']);
fNameCell=extractfield(files,'name');
files = files(~contains(fNameCell,'config') & ~contains(fNameCell,'min.mat'));
Nf    = length(files);
%
% initialize counter...
fprintf('\n \n')
ensNum = 1;
waves  = struct();
for ii= 1:Nf
    %
    % load the raw data
    inFileName = sprintf([filePrefix,'%03d.mat'],ii);
    fin  = [outDir,inFileName];
    fprintf(['loading file:   %s \n'],inFileName)
% $$$     in = load(fin,'VelEast','VelNorth','VelUp1','VelUp2','VelBeam5','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','mab');
    in = load(fin,'Velocity_East','Velocity_North','Velocity_Up','qcFlag','HeadingOffset','Time','Heading','Pitch','Roll','Pressure','bin_mab');
    %
    if ii==1
        mab = in.bin_mab;
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
        QC   = in.qcFlag(:,(jj-1)*Na+1:Na*jj);
        fGood= sum(QC,2)/Na;
        % only use bins in water >99 % of the time
        bins = find(fGood>0.99);
        Nb   = length(bins);
        %
        % mean time to nearest minute
        t = in.Time((jj-1)*Na+1:Na*jj);
        tavg = round(mean(t)*1440)/1440;
        t = (t-t(1))*86400;
        %
        waves.Time(1,ensNum) = tavg;
        if Nb<2
            disp('not enough good data')
            waves.Hs(1,ensNum) = nan;
            waves.Tm(1,ensNum) = nan;
            waves.Suu(:,ensNum)= nan;
            waves.Svv(:,ensNum)= nan;
            waves.Spp(:,ensNum)= nan;
            waves.Spu(:,ensNum)= nan;
            waves.Spp(:,ensNum)= nan;
            waves.Spv(:,ensNum)= nan;
            waves.mean_dir(1:2,ensNum)   =nan;
            waves.mean_spread(1:2,ensNum)=nan;
            waves.mSxx(1,ensNum) = nan;
            waves.mSxy(1,ensNum) = nan;
            waves.mSyy(1,ensNum) = nan;
            waves.Z2(:,ensNum) = nan;
            ensNum = ensNum+1;
            continue
        end
        %
        % allocate current variables: (note, bad QC are zeros not nans, and not interpolated)
        U = in.Velocity_East (bins,(jj-1)*Na+1:Na*jj);
        V = in.Velocity_North(bins,(jj-1)*Na+1:Na*jj);
        W = in.Velocity_Up   (bins,(jj-1)*Na+1:Na*jj);
        P = in.Pressure(2,(jj-1)*Na+1:Na*jj)*1e4/(rho0*g);
        %
        % depth of pressure sensor
        dpthP= mean(P,2);% not accounting for sensor height above bed
        % assuming sensor is at "hab" above bed level        
        dpth = dpthP + hab;
        % depth of bins
        dpthU= dpthP-mab(bins);
        %
% $$$         % Use WAFOS package:
% $$$         % first pass, just use (U,V,W,P):
% $$$         xn  = [t', (P'-dpthP)*1e4, U', V', W'];
% $$$         pos0= [0, 0,-dpthP, 9, 1];
% $$$         posU= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*10, ones(Nb,1)];
% $$$         posV= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*11, ones(Nb,1)];
% $$$         posW= [zeros(Nb,1),zeros(Nb,1),-dpthU,ones(Nb,1)*12, zeros(Nb,1)];        
% $$$         pos = [pos0;posU;posV;posW];
% $$$         [Sd,D,Sw,Fcof,Gwt,Sxy,Sxy1] = dat2dspec_DG(xn,pos,dpth,Ne/4,90,'MLM');% 256 frequency bins & 90 directional bins
% $$$         %
% $$$         % quick plot in polar coords
% $$$         x = Sd.f'.*cos(Sd.theta);
% $$$         y = Sd.f'.*sin(Sd.theta);
% $$$         figure, contourf(x,y,log10(Sd.S),[-4:0.05:3],'edgecolor','none'),colormap(cmocean('thermal'))
% $$$         %
% $$$         %
% $$$         % estimate mean direction
% $$$         df = Sd.f(2)-Sd.f(1);
% $$$         dT = D.theta(2)-D.theta(1);
% $$$         I  = find(Sd.f>=1/20 & Sd.f<=1/4);
% $$$         A1 = sum( D.S.*cos(D.theta)*dT, 1);
% $$$         B1 = sum( D.S.*sin(D.theta)*dT, 1);
% $$$         %
        %
        % estimate bulk wave statistics on ~7 min chuncks, 2/3 overlap
        [Suu,fq]=welch_method(U',dt,chnks,olap);
        [Svv,fq]=welch_method(V',dt,chnks,olap);
        [Spp,fq]=welch_method(P',dt,chnks,olap);
        [Suv,fq]=welch_cospec(U',V',dt,chnks,olap);        
        [Spu,fq]=welch_cospec(repmat(P',1,Nb),U',dt,chnks,olap);
        [Spv,fq]=welch_cospec(repmat(P',1,Nb),V',dt,chnks,olap);
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
        %
        % convert pressure to surface elevation
        cP     =  cosh(k.*dpth)./ cosh(k.*(dpth-dpthP)); % SePP=Spp.*cP.^2
        % convert velocity to surface elevation
        cU2eta = (om./(g*k)).*cosh(k.*dpth)./cosh(k.*(dpth-dpthU'));% SeUU = SeUU.*cU.^2        
        % convert velocity to surface velocity
        cU2srf =              cosh(k.*dpth)./cosh(k.*(dpth-dpthU'));% SeUU = SeUU.*cU.^2        
        %
        % estimate max resolvable frequency where depth bin is below wave-length/(pi)
        fmax_vs_bin = fq(Ne/2-sum(dpthU'>=1./(k)))';
        % nan surface transformed fields above fmax
        %        Spp(fq>fmax(1))=nan;
        for kk=1:Nb
        Suu(fq>fmax_vs_bin(kk),kk)=nan;
        Svv(fq>fmax_vs_bin(kk),kk)=nan;
        Suv(fq>fmax_vs_bin(kk),kk)=nan;
        Spu(fq>fmax_vs_bin(kk),kk)=nan;
        Spv(fq>fmax_vs_bin(kk),kk)=nan;
        end
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
        % estimate best guess for surface elevation spectrum
% $$$         % first, average the velocity spectra
% $$$         SeAVG = nanmean(SeUU+SeVV,2);
% $$$         % use z-score as a weighting function
% $$$         norm   = [tanh(pi*SePP./SeAVG), ones(Ne/2,1)];
% $$$         norm(2:end-1,:) = (norm(1:end-2,:) + norm(2:end-1,:) + norm(3:end,:))/3;
% $$$         SeBest = nansum([SeAVG SePP].*norm,2)./sum(norm,2);
        %
        if hour(tavg)==12
            fig1 = figure;
            plt = semilogy(fq,nanmean(SeUU+SeVV,2),'-r',fq,SePP,'--k','linewidth',3);
            set(gca,'xlim',[1/(30*60) 1/3])
            grid on
            xlabel('cycles per second','interpreter','latex')
            ylabel('m$^2$/Hz','interpreter','latex')
            title(datestr(tavg),'interpreter','latex')
            h = legend(plt,'$S_{uu}+S_{vv}$','$S_{\eta\eta}$');
            set(h,'location','northeast','interpreter','latex')
            figname = sprintf('%s/figures/%s_spectra_%s.pdf',outRoot,L1FRoot,datestr(tavg,'mmddHHMM'));
            if ~exist([outRoot,filesep,'figures'])
                eval(['!mkdir -p ',[outRoot,filesep,'figures/']])
            end
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
        I = find(fq>=1/20 & fq<=1/4);
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
            waves.frequency  = fq(1:I(end))';
            waves.wavenumber = k(1:I(end))';            
        end
        waves.Hs     (1,ensNum)   = Hs;
        waves.Tm     (1,ensNum)   = Tm;
        waves.depth  (1,ensNum)   = dpth;
        waves.Spp    (:,ensNum)   = SePP(1:I(end));
        waves.Suu    (:,ensNum)   = nanmean(SUU(1:I(end),:),2);
        waves.Svv    (:,ensNum)   = nanmean(SVV(1:I(end),:),2);
        waves.Spu    (:,ensNum)   = nanmean(SPU(1:I(end),:),2);
        waves.Spv    (:,ensNum)   = nanmean(SPV(1:I(end),:),2);
        waves.mdir   (1:2,ensNum) = [nanmean(mdir1)   ,nanmean(mdir2)];
        waves.mspread(1:2,ensNum) = [nanmean(mspread1),nanmean(mspread2)];
        waves.a1     (:,ensNum)   = nanmean(a1,2);
        waves.b1     (:,ensNum)   = nanmean(b1,2);        
        waves.a2     (:,ensNum)   = nanmean(a2,2);
        waves.b2     (:,ensNum)   = nanmean(b2,2);        
        waves.mSxx   (1,ensNum)     = nanmean(mSxx);
        waves.mSxy   (1,ensNum)     = nanmean(mSxy);
        waves.mSyy   (1,ensNum)     = nanmean(mSyy);
        waves.Z2     (:,ensNum)    = Z2;
        ensNum = ensNum+1;
    end
    %
end
fprintf(['saving: %s \n'],[L1dir,L1FRoot])
if ~exist(L1dir,'dir')
    eval(['!mkdir -p ',L1dir])
end
save([L1dir,L1FRoot,'.mat'],'currents','waves')
fprintf('done! \n')