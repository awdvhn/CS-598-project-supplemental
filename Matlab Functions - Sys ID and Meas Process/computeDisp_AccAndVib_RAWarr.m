function [time, F_deNoise, x_acc, v_acc, a_acc, x_vib, v_vib, d, c,...
    iFstart, iFend] = computeDisp_AccAndVib_RAWarr(timeRAW, Fraw_ham, aRAW_acc, vRAW_vib )

% Denoise Hammer Force, Get Impulse Start and End instants after denoise
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(Fraw_ham); 

% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

% Trim Time Instants to after Hit Region
time = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

% Trim Denoised Force to after Hit Region
F_deNoise = F_deNoise(I_ForceStart);

% Trim Vibrometer data to after Hit Region
v_vib = detrend(vRAW_vib(I_ForceStart));

% Trim Accelerometer data to after Hit Region
a_acc = detrend(aRAW_acc(I_ForceStart));

%% Create Filter for Processing
Fs = 1/mean(diff(time));
Fnyq = Fs/2;

% Approximate of natural frequency
[f_FFTacc, aFFT_acc] = regFFT(time, a_acc); 
[~, ires_acc] = max(aFFT_acc(f_FFTacc > 5));
fres_acc = f_FFTacc(ires_acc);

% L-Pass Filter for Experimental Noise
fc_LowPass = 2*fres_acc;
Fc_LowPass = fc_LowPass/Fnyq;
[b,a] = butter(5,Fc_LowPass,'low');

% High-Pass Filter for Integration
fc = fres_acc*0.35;
Fc = fc/Fnyq; % Filter cutoff ratio
[d,c] = butter(3,Fc,'high');

%% Integrate to get displacement

% Vibrometer Data -> Displacement
v_vib = filtfilt(b, a, v_vib);

x_vib = detrend(cumtrapz(time, v_vib));
x_vib = filtfilt(d, c, x_vib);

% Accelerometer Data -> Displacement 
a_acc = filtfilt(b, a, a_acc);

v_acc = detrend(cumtrapz(time, a_acc));
v_acc = filtfilt(d, c, v_acc);

x_acc = detrend(cumtrapz(time, v_acc));
x_acc = filtfilt(d, c, x_acc);

end

