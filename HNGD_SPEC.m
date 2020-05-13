function [HNGD,HNGD_stat]=HNGD_SPEC(wpt)

%addpath('/home/user/Desktop/Speech_codes/Speech_codes/GRPDEALY/voicebox')
% [inputsignal, fs]=audioread(wavname);
% inputsignal=resample(inputsignal,8000,fs);fs=8000;

[inputsignal,fs]=audioread(wpt);
f_len=30;f_shft=15;w2=hamming(f_len*fs/1000);

eps=10^-12;
mfccsignal=inputsignal+eps;
mfccsignal = filter([1 -0.97], 1, mfccsignal);

FrameLenMillisec    = f_len;
FrameShiftMillisec  = f_shft;

FrameLen = round((FrameLenMillisec/1000)*fs);
NFFT            = 512;
FrameShift = round((FrameShiftMillisec/1000)*fs);
Frames = enframe(mfccsignal, w2', FrameShift); %%% NOTE: NO WINDOWING
WindowedFrames=(Frames .* repmat(hamming(size(Frames,2))', size(Frames,1), 1))';
% WindowedFrames=Frames;

frames=WindowedFrames';
[frame_num, frame_length] = size(frames);

delay_vector = [1:1:frame_length];
delay_matrix = repmat(delay_vector, frame_num, 1);

delay_frames = frames .* delay_matrix;

x_spec = fft(frames', NFFT);
y_spec = fft(delay_frames', NFFT);
x_spec = x_spec(1:NFFT/2, :);
y_spec = y_spec(1:NFFT/2, :);

NGD=(real(x_spec).*real(y_spec) + imag(y_spec) .* imag(x_spec));
HNGD=abs(hilbert(NGD));

num_coeff=20;
grp_phase=NGD;
grp_phase(isnan(grp_phase)) = 0.0;

cep = dct(grp_phase);
cep = cep(2:num_coeff, :)';
f1=mean(HNGD');
                f2=std(HNGD');
                f3=skewness(HNGD');
                f4=kurtosis(HNGD');
                f5=range(HNGD');
                f6=mean(abs(diff(HNGD')));
                HNGD_stat=[f1 f2 f3 f4 f5 f6];

end