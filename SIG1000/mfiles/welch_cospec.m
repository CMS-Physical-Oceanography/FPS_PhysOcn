function [CoSP,fr,QuSP,COH,PHI]=welch_cospec(data0,data1,dt,M,lap);
%
% USAGE: [CoSP,fr,QuSP,fr,Coh,Phi]=welch_cospec(data0,data1,dt,M,lap);
%
% function estimates CoSPectral density of "data0" and "data1" using "M" subwindows,
% with "lap" fractional overlap. The dof = 2*M for lap=0,

% size of input data
sd = size(data0); % time must run down
% get the number of samples in each chunk
Ns=floor(sd(1)/(M-(M-1)*lap));
M =floor(M);
% get how many indices to change by in the loop
ds=floor(Ns*(1-lap));
ii=1:Ns;
%
% use hanning window
win=hanning(Ns);
win=repmat(win,[1 sd(2)]);
%
if iseven(Ns)
    inyq = 1;% don't double the nyquist frequency
    stop = Ns/2+1;
else
    inyq = 0;%
    stop = (Ns+1)/2;
end
%
fr = [0:stop-1]'/(dt*Ns);
%
SXX=zeros(stop,sd(2));
SYY=zeros(stop,sd(2));
CXY=zeros(stop,sd(2));
QXY=zeros(stop,sd(2));
for m=1:M
    inds   = ii+(m-1)*ds;% indices in current block
    x  = data0(inds,:);% data0 from ...
    y  = data1(inds,:);% data0 from ...    
    sx2i= var(x);% input variance,
    sy2i= var(y);% input variance,     
    x  = win.*detrend(x);% detrend and window
    y  = win.*detrend(y);% detrend and window        
    sx2f= var(x);% reduced variance
    sy2f= var(y);% reduced variance    
    x  = repmat((sx2i./sx2f).^.5,[Ns 1]).*x;
    y  = repmat((sy2i./sy2f).^.5,[Ns 1]).*y;    
    %
    X  = fft(x); X  = X(1:stop,:);% positive frequencies
    Y  = fft(y); Y  = Y(1:stop,:);% positive frequencies
    %
    AXX = X.*conj(X);% amplitude squared
    AXX(2:end-inyq,:)  = 2*AXX(2:end-inyq,:);
    AYY = Y.*conj(Y);% amplitude squared
    AYY(2:end-inyq,:)  = 2*AYY(2:end-inyq,:);
    AXY = X.*conj(Y);
    AXY(2:end-inyq,:)   = 2*AXY(2:end-inyq,:);
    %
    SXX  =SXX+AXX;
    SYY  =SYY+AYY;
    CXY  =CXY+AXY;
end
SXX=SXX*dt/M/Ns;
SYY=SYY*dt/M/Ns;
CXY=CXY*dt/M/Ns;

CoSP= real(CXY);
QuSP= imag(CXY);
COH = abs(CXY)./sqrt(SXX.*SYY );
PHI = atan2(-QuSP,CoSP);


