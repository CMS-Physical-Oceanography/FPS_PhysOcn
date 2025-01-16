function [psd,fr]=welch_method(data,dt,M,lap);
%
% USAGE: [psd,fr]=welch_method(data,dt,M,lap);
%
% function estimates power spectral density of "data" using "M" subwindows,
% with "lap" fractional overlap. The dof = 2*M for lap=0,
% size of input data
sd = size(data); % time must run down
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
SX=zeros(stop,sd(2));
for m=1:M
    inds   = ii+(m-1)*ds;% indices in current block
    x  = data(inds,:);% data from ...
    s2i= var(x);% input variance, 
    x  = win.*detrend(x);% detrend and window    
    s2f= var(x);% reduced variance
    x  = repmat((s2i./s2f).^.5,[Ns 1]).*x;
    %
    X  = fft(x);
    X  = X(1:stop,:);% positive frequencies
    A2 = X.*conj(X);% amplitude squared
    A2(2:end-inyq,:)  = 2*A2(2:end-inyq,:);
    SX=SX+A2;
end
psd=SX*dt/M/Ns;
