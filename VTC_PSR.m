function [VTC, PSR, PSR_stat] = VTC_PSR(wnm)
clc;

% addpath('/home/user/Desktop/Speech_codes/Speech_codes/GRPDEALY/voicebox');
%% have to noramlize each and every file

nlpc=10;
VTC=[];
PSR=[];

        [wa,fs]=audioread(wnm);
        wa=resample(wa,8000,fs);fs=8000;
        wa=wa/(abs(max(wa)));
        wa=smooth(wa);
        [~,~,~,zfe,~,GCI,VAD]=EPOCH_SOE_F0_ZFF(wa,fs);

        wa=wa.*VAD;
        zfe_n=zfe.*VAD;
        xy   = dot(abs(wa),abs(zfe));
nx   = norm(wa);
ny   = norm(zfe);
nxny = nx*ny;
Cs   = xy/nxny;
VTC=Cs;
%         
        %% 50 ms frames lpcs

        le=20*fs/1000; shft=10*fs/1000;
        [L,~]=LPC_RESIDUAL(wa,le,shft,nlpc);
        %       
        sl=length(L);
        ng=0;
        for gnum=1:length(GCI)
            
            if(VAD(GCI(gnum))==1)
            s_win=L(max(1,GCI(gnum)-40):min(sl,GCI(gnum)+40));
            [val,loc]=max(abs(s_win));
            n_win= L(max(1,loc-16):min(sl,(loc+16)));
            [pks,locs]=findpeaks(n_win);
            mn=(sum(abs(pks))-val)/(length(pks)-1);
            ng=ng+1;
            PSR(ng)=val/mn;
            end
        end


%%
T=PSR';
 f1=mean(T);
                f2=std(T);
                f3=skewness(T');
                f4=kurtosis(T');
                f5=range(T);
                f6=mean(abs(diff(T)));
                FT_stat=[f1 f2 f3 f4 f5 f6];
                PSR_stat=FT_stat;
end
