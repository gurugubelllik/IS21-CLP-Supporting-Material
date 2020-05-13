function [S,cep] = LPCC(inputsignal, fs,p,alpha)

% [inputsignal, fs]=audioread(wavname);
% inputsignal=resample(inputsignal,8000,fs);fs=8000;

eps=10^-12;
mfccsignal=inputsignal+eps;
mfccsignal = filter([1 -0.97], 1, mfccsignal);

FrameLenMillisec    = 50;
FrameShiftMillisec  = 10;
num_coeff=20;

FrameLen = round((FrameLenMillisec/1000)*fs);
nfft            = 512;
FrameShift = round((FrameShiftMillisec/1000)*fs);
Frames = enframe(mfccsignal, hann(FrameLen), FrameShift); %%% NOTE: NO WINDOWING
WindowedFrames=(Frames .* repmat(hamming(size(Frames,2))', size(Frames,1), 1))';
frames=WindowedFrames;

X = fft(frames,nfft);
R = ifft((abs(X)).^(2.*alpha)); % LP-alpha paper my TASLP-2016 paper
R = R./size(frames,1); % Biased autocorrelation estimate
countStablized=0;
for i = 1:size(R,2)
    r = R(2:p+1,i);
    Autocorr = toeplitz(R(1:p,i));
    a2 = inv(Autocorr)*r;
    a(:,i) = [1;-a2];
    
    % LP inverse filter spectra
    ifspec(:,i) = abs(fft(a(:,i),nfft));
end

S = 1./(ifspec.^2);

cep = dct(log(S));
cep = cep(2:num_coeff+1, :)';