clc
close all
clear all
wavpath='aaa.wav';


%% Baseline Specral and prosody feature extraction

cqcc_features=CQCC_Feat(wavpath); %% outs N X D features, N: no. of frames D: 39 dimension


MFCC_features=MFCC_Feat(wavpath); %% outs N X D features, N: no. of frames D: 39 dimension


SFFB_LTAS_features=LTAS_SFF_Feat(wavpath); %% outs 1 X 2019 feature VECTOR

[Zff,Lpr,Gvv,Zff_st,Lpr_st,Gvv_st]=extract_spcom_excitaion_feats_frm_given_wpt(wavpath); % Zff,LPR,GVV : number of frames X 12.; STAT 1X(12*6) 

[HNGD,HNGD_stat]=HNGD_SPEC(wavpath); % HNGD 256 X N : no of frames, HNGD_stat 1X(256*6)

[VTC, PSR, PSR_stat] = VTC_PSR(wavpath); % VTC 1X1 , PSR 1 X number of epoch locations
