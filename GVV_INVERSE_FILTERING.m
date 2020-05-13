function [gvv]=GVV_INVERSE_FILTERING(wav,fs,vnv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used for find Glotal volume velocity wave from speech Signal
% Edited by: G Krishna  (Speech Lab, IIIT Hyderabad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------       inputs         --------
% wav       : input raw speech signal
% fs        : sampling freqnecy
% vnv       : voiced unvoiced region
% --------      outputs         --------
% gvv       : glottal volume velocity found by inverse filtering method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wav=resample(wav,8000,fs);

fs=8000;

if nargin<3
    [epoch,soe,F0,zfSig,vgci,gci,vad]=EPOCH_SOE_F0_ZFF(wav,fs);    
    vnv=vnvseg(wav,zfSig,fs,10);    
end

analysis_window_length=20*fs/1000; % 20msec or 30msec
wav=wav./max(abs(wav));
wav_seg=buffer(wav,analysis_window_length);

% wind=hamming(analysis_window_length);
% for i=1:size(wav_seg,2)
% wav_seg(:,i)=wav_seg(:,i).*wind;
% end

[a,r]=lpc(wav_seg,10);

for i=1:size(wav_seg,2)
    lp_res(:,i)=filter(a(i,:),1,wav_seg(:,i));
end

lp_res=reshape(lp_res,size(lp_res,1)*size(lp_res,2),1);

lp_res=smooth(lp_res,round(fs/1000));

lp_res=cumsum(lp_res);

gvv=remTrend(lp_res,round(10*fs/1000));

gvv=gvv./max(abs(gvv));

gvv=gvv(1:length(wav)).*vnv;

% figure
% a=subplot(211);
% plot(gvv);hold on;plot(vnv,'r');
% b=subplot(212);
% plot(wav)
% linkaxes([a,b],'x')
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

[maxv, maxi]=max(out(minPitch*fs/1000:maxPitch*fs/1000,:));

x=(minPitch:0.5:maxPitch)*fs/1000+2;
pLoc=maxi(vad>frameSize*0.8)+minPitch*fs/1000;
y=hist(pLoc,x);
y=y/length(pLoc);

[val, idx]=max(y);
idx=round(idx/2)+minPitch;
end


function [out]=remTrend(sig,winSize)

window=ones(winSize,1);
rm=conv(sig,window);
rm=rm(winSize/2:length(rm)-winSize/2);

norm=conv(ones(length(sig),1),window);
norm=norm(winSize/2:length(norm)-winSize/2);

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


function [vnvsig] = vnvseg(s,zf,fs,nmean)

if(~exist('nmean'))
    nmean=10;
end;

s=s(:);
ds=diff(s);
ds(end+1)=ds(end);

s=s/sqrt(sum(s.^2)/length(s));
ds=ds/sqrt(sum(ds.^2)/length(ds));
zf=zf/sqrt(sum(zf.^2)/length(zf));
%
se=RunMean(ds.^2,nmean*fs/1000);
zfe=RunMean(zf.^2,nmean*fs/1000);

zfbys=zfe./se;
zfbysevi=-tansig(zfbys-10);
zfevi=1-exp(-10*zfe);

vnvfevi=zfevi;
%vnvevi = [zfevi(:) zfbysevi(:) rse10evi(:)];
vnvsig=vnvfevi > 0.5;
vnvsig=double(vnvsig);
return;

end



function [yavg, ysum] = RunMean(sig, N, wintype)

if( ~ exist('wintype'))
    wintype	= 'RECT';
end

M	= length(sig);
Nby2	= floor(N/2);
Ntail	= N - Nby2 - 1;

switch wintype
    case {'REC', 'RECT', 'rec', 'rect', '1'}
        h	= ones(N,1);
        
    case {'HAM', 'HAMM', 'ham', 'hamm', '2'}
        h	= hamming(N);
        
    case {'HAN', 'HANN', 'han', 'hann', '3'}
        h	= hanning(N);
        
    otherwise
        disp('Error : Unknown window type!!!');
        exit(1);
end

x	= conv(sig,h);
ysum	= x(Nby2+1:M+N-1-Ntail);

xdiv	= conv(ones(size(sig)),h);
x	= x ./ xdiv;
yavg	= x(Nby2+1:M+N-1-Ntail);

return;

end