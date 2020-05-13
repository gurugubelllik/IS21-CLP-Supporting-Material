function  [Zff,Lpr,Gvv,Zff_st,Lpr_st,Gvv_st]=extract_spcom_excitaion_feats_frm_given_wpt(wpt)
% smt_y 0 no smooth, 1 smooth , 2 smooth after vad
% aaa=bbb;
%mn_fld='/home/user/Desktop/Trails/data_preparation/final_data_all3/ALL_T';
%addpath('/home/user/Desktop/Speech_codes/Speech_codes/GRPDEALY/voicebox');
%addpath('/home/user/Desktop/Trails/data_preparation/final_data_all3/set-1');

%% have to noramlize each and every file

smt_y=1;vad_y=1;

FT_zff=[];FT_lp=[];FT_gvv=[];


uttr_lbls=[];
fls_pr_spkr=[];
fls_pr_file=[];
r=0.98; %% zff
nlpc=10;%% lp coeffs


 [wa,fs]=audioread(wpt);
wa=resample(wa,8000,fs);fs=8000;
wa=wa/(abs(max(wa)));
if smt_y==1
    wa=smooth(wa);
end
[~,~,~,~,~,~,VAD]=EPOCH_SOE_F0_ZFF(wa,fs);

if vad_y==1
    wa=wa.*VAD;
end
if smt_y==2
    wa=smooth(wa);
end

G=GVV_INVERSE_FILTERING(wa,fs);
[WinLen]=xcorrWinLen(wa,fs);

[Z,zp_zff]=zeroFreqFilter_1(wa,fs,WinLen,r);
lpcs=LPCC_Feat(wa,fs,20,1);
le=20*fs/1000; shft=10*fs/1000;
[L,~]=LPC_RESIDUAL(wa,le,shft,nlpc);
L=abs(L);
%%
% G=filter([1,0.98],1,L);
%%
GM=melcepst(G,fs,30);reshape(GM',1,[]);
ZM=melcepst(Z,fs,30);
LM=melcepst(L,fs,30);
Zff=ZM;Lpr=LM;Gvv=GM;
Zff_st=FT_stat(Zff);Lpr_st=FT_stat(Lpr);Gvv_st=FT_stat(Gvv);
end
function FT_st=FT_stat(T)
f1=mean(T);
f2=std(T);
f3=skewness(T);
f4=kurtosis(T);
f5=range(T);
f6=mean(abs(diff(T)));
FT_st=[f1 f2 f3 f4 f5 f6];
end
