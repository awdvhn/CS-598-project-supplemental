function [f,ampFFT, phaseFFT] = regFFT(Time,y)
% Summary of this function goes here
%   Detailed explanation goes here

L = length(y); % length of hanning window = length of data
nfft = 2^nextpow2(L); % Transform length
y_HannWnd = y; % No hanning window applied
Ydft_HannWnd = fft(y_HannWnd,nfft)/L;

% at all frequencies except zero and the Nyquist
mYdft = abs(Ydft_HannWnd);
mYdft = mYdft (1:nfft/2+1);
mYdft (2:end-1) = 2* mYdft(2:end-1); % symmetric fft factor
dt = diff(Time);
Fs = 1/dt(1);
f = Fs/2*linspace(0,1,nfft/2+1); 
ampFFT = mYdft; %no window correction factor

% Compute angle
Ydft_HannWnd = fftshift(Ydft_HannWnd);
phaseFFT = angle(Ydft_HannWnd);
phaseFFT = phaseFFT(nfft/2:end);

end



