function [freqAVG_Hz, zetaAVG, Xmax, freq_tFILT, zeta, logXratio_t,...
    xA, phase_t, freq_t, X_tFILT, X_t, b, a, iXmax] = freeLinSDoF_SysID_Hilbert(tspan, x)
% Algorithm:
%           1) Computes Analutical Solution by the mean of Hilbert
%           transform called xA
%           2) Compute instantanneous phase, frequency, envelope,
%           normalized envelope log decrement using xA
%           3) Create low-pass using 1/10 x average frequency of HT
%           4) apply filter instantanneous frequency and normalized log
%           decrement
%           5) Compute average frequency
%           6) Compute zeta(t) using avergae frequency and then the
%           average of zeta
%
%           NOTE: averaging is done on halve the points at the center of
%           duration considered i.e. 
%                   rangeOi = length(X_t)/4):ceil(3*length(X_t)/4

%% Compute analytical solution of given displaceement by Hilbert Transform
xA = hilbert(x);

%% Compute the Frequency
dt = mean(diff(tspan));
Fs = 1/dt;
phase_t = unwrap(angle(xA)); % [rad]
freq_t = Fs*gradient(phase_t)/(2*pi); % [Hz]

%% Compute Envelope amplitude
X_t = abs(xA);

%% Create a low-pass filter with 1/10 the approximate oscillations of interest 
rangeOi = floor(length(X_t)/3):ceil(2*length(X_t)/3);
fnApp = mean(freq_t(rangeOi));
Fc = fnApp/10/(Fs/2);
[b,a] = butter(1,Fc,'low');

%% Filter Frequency and Evelop Log Decrement
freq_tFILT = filtfilt(b,a,freq_t);
X_tFILT = filtfilt(b,a,X_t);

%% Compute Log Decrement
[Xmax, iXmax] = max(X_tFILT);
logXratio_t = log(X_tFILT./Xmax);

%% Compute average frequency, inst zetas and avergae zeta
freqAVG_Hz = mean(freq_tFILT(rangeOi));
zeta = - logXratio_t./(2*pi*freqAVG_Hz * (tspan - tspan(iXmax)));
zetaAVG = mean(zeta(rangeOi));

end

