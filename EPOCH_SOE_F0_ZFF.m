function [epoch,sgci,F0,zfSig,vgci,gci,vad]=EPOCH_SOE_F0_ZFF(wav,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used for finding strengths of excitation from ZFR Signal
% Authors: B. yegnanarayana and K. S. R. Murthy
% Edited by: G Krishna  (Speech Lab, IIIT Hyderabad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------       inputs         --------
% wav       : input raw speech signal
% fs        : sampling freqnecy
% --------      outputs         --------
% epoch     : epoch signal in voiced region for a speech signal
% soe       : stremgth of excitation os speech signal
% F0        : pitch/F0 contour of a speech signal
% zfSig     : zero frequency filtered signal
% vgci      : glottal closer in vloced regions
% gci       : estimated glottal closer instants for complete signal
% vad       : voice activity region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [wav,fs]=wavread(filename);
wav=wav./max(abs(wav));
winLength=xcorrWinLen(wav,fs);
% % 
winLength=winLength;
% winLength=10;
zfSig=zeroFreqFilter(wav,fs,winLength);

vad=getVad(wav,fs);
vad=durational_constraint(vad,30*fs/1000);

gci=epoch_extract(wav,zfSig);

t=zeros(length(wav),1);
t(gci)=1;
vgci=find(vad.*t); %vgci is gci's in voiced regions only

epoch=zeros(size(wav));
epoch(vgci)=1;

sgci=abs(zfSig(gci)-zfSig(gci-1));
sgci=sgci/max(abs(sgci));


fgci=diff(gci);
fgci(end+1)=fgci(end);
F0=fs./fgci;

% close all;
% figure
% a=subplot(311);
% plot(wav./max(abs(wav)));grid
% b=subplot(312);
% plot(zfSig./max(abs(zfSig)));hold on;plot(vad,'r');plot(soe*-0.9,'r');grid
% c=subplot(313);
% temp1=(zfSig./max(abs(zfSig)));
% temp2=wav./max(abs(wav));
% temp3=zeros(size(wav));
% temp3(temp2~=0)=temp1((temp2~=0))./temp2(temp2~=0);
% 
% plot(smooth(temp3,20*fs/1000));
% linkaxes([a b c],'x')
% grid
% 
% figure
% a=subplot(311);
% plot(wav./max(abs(wav)));grid on;
% b=subplot(312);
% plot(zfSig./max(abs(zfSig)));hold on;plot(vad,'r');plot(soe*-0.9,'r');grid on;
% c=subplot(313);
% plot(F0,'*');grid on;
% linkaxes([a b c],'x')
end


function [idx]=xcorrWinLen(wav,fs)

frameSize=30*fs/1000;
frameShift=20*fs/1000;

en=conv(wav.^2,ones(frameSize,1));
en=en(frameSize/2:end-frameSize/2);
en=en/frameSize;
en=sqrt(en);
en=en>max(en)/5;

b=buffer(wav,frameSize,frameShift,'nodelay');
vad=sum(buffer(en,frameSize,frameShift,'nodelay'));

FUN=@(x) xcorr((x-mean(x)).*hamming(length(x)),'coeff')./xcorr(hamming(length(x)),'coeff');
out=blkproc(b,[frameSize,1],FUN);

out=out(frameSize:end,:);

minPitch=3;  %2 ms == 500 Hz.
maxPitch=16; %16 ms == 66.66 Hz.

[~, maxi]=max(out(minPitch*fs/1000:maxPitch*fs/1000,:));

x=(minPitch:0.5:maxPitch)*fs/1000+2;
pLoc=maxi(vad>frameSize*0.8)+minPitch*fs/1000;
y=hist(pLoc,x);
y=y/length(pLoc);

[~, idx]=max(y);
idx=round(idx/2)+minPitch;
end


function [zfSig]=zeroFreqFilter(wav,fs,winLength)
% Difference the speech signal...
dwav=diff(wav);
dwav(end+1)=dwav(end);
dwav=dwav/max(abs(dwav))/1.001;
N=length(dwav);

% Pass the differenced speech signal twice through zero-frequency resonator..

zfSig=cumsum(cumsum(cumsum(cumsum(dwav))));

% Remove the DC offset introduced by zero-frquency filtering..
winLength=round(winLength*fs/1000);
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig(N-winLength*3:N)=0;
end

function [out]=remTrend(sig,winSize)

window=ones(winSize,1);
rm=conv(sig,window);
rm=rm(winSize/2:length(rm)-winSize/2);

norm=conv(ones(length(sig),1),window);
norm=norm(winSize/2:length(norm)-winSize/2);
% yyy=size(norm)
zzz=size(rm);
rm=rm./norm;
out=sig-rm;
end

function [vad]=getVad(sig,fs)

winLength=20*fs/1000;
en=conv(abs(sig),ones(winLength,1));
en=en(winLength/2:length(en)-winLength/2);
en=en/max(en);

%figure; plot(sig); hold on; plot(en/max(abs(en)),'r');

vad=en>0.1;
end


function[in_reg]=durational_constraint(in_reg,dur_cons)

% dur_cons = durration in samples

in_reg(1)=0;in_reg(end)=0;
dp=diff(in_reg);
pos=find(dp>0);
neg=find(dp<0);
if length(pos)~=length(neg)
    if(neg(1)<pos(1))
        neg=neg(2:end);
    end
    if(neg(end)<pos(end))
        pos=pos(1:end-1);
    end
end

if length(pos)==length(neg)
    for i=1:length(pos)
        dur=length(pos(i):neg(i));
        if (dur<dur_cons)
            in_reg(pos(i):neg(i))=0;
        end
        
    end
end
end


function [gci]=epoch_extract(wav,zfSig)

[pzc, pslope]=zerocros(zfSig,'p');
[nzc, nslope]=zerocros(zfSig,'n');

l=min(length(pzc),length(nzc));

dslope=abs(pslope(1:l))-abs(nslope(1:l));
if(sum(dslope>0)>length(dslope)/2)
    gci=pzc;
    %		disp('Polarity of the signal: +ve');
else
    gci=nzc;
    wav=-wav;
%     zfSig=-zfSig;
    %		disp('Polarity of the signal : -ve');
end;

%	p=zeros(length(wav)-1,1); p(gci)=1; p=p.*vad; gci=find(p==1); %uncomment if u want instants in voiced regions only
p=zeros(length(wav)-1,1); p(gci)=1; gci=find(p==1);

gci=gci+3;

end


function [f,s]=zerocros(x,m)

if nargin<2
    m='b';
end

s=x>=0;
k=s(2:end)-s(1:end-1);

if any(m=='p')
    f=find(k>0);
elseif any(m=='n')
    f=find(k<0);
else
    f=find(k~=0);
end

s=x(f+1)-x(f);

end