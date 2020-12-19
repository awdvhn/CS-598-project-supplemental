function [x_acc, v_acc, a_acc, t_acc, x_vib, v_vib, t_vib, d, c, t_ham, F_deNoise, aRAW_acc, vRAW_vib, iFstart, iFend] = computeDisp_AccAndVib(expDataFile)
%% Load Data
load(expDataFile)

% Hammer Force: Get Maximum and instant after denoise
t_ham = MAT_006_Time_record_Hammer_(:,1);
F_ham = MAT_006_Time_record_Hammer_(:,2);
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(F_ham); % Eliminate useless measurements
                                                              % of the hammer
I_ForceStart = iFstart:length(t_ham);

% Vibrometer data
t_vib = MAT_008_Time_record_Vibro__(:,1); t_vib = t_vib(I_ForceStart) - t_vib(I_ForceStart(1));
vRAW_vib = detrend(MAT_008_Time_record_Vibro__(:,2));
v_vib = vRAW_vib(I_ForceStart);

% Accelerometer data
t_acc = MAT_007_Time_record_Acc__(:,1); t_acc = t_acc(I_ForceStart) - t_acc(I_ForceStart(1));
aRAW_acc = detrend(MAT_007_Time_record_Acc__(:,2)); 
a_acc = aRAW_acc(I_ForceStart);

%% Integrate to get displacement

% Create High-Pass Filter
Fs = 1/mean(diff(t_acc));
Fnyq = Fs/2;

[f_FFTacc, aFFT_acc] = regFFT(t_acc, a_acc); % Approximate of natural frequency
[~, ires_acc] = max(aFFT_acc(f_FFTacc > 5));
fres_acc = f_FFTacc(ires_acc);
fc = fres_acc*0.35;

Fc = fc/Fnyq; % Filter cutoff ratio
[d,c] = butter(3,Fc,'high');

% Vibrometer Data -> Displacement
x_vib = detrend(cumtrapz(t_vib, v_vib));
x_vib = filtfilt(d, c, x_vib);

% Accelerometer Data -> Displacement
v_acc = detrend(cumtrapz(t_acc, a_acc));
v_acc = filtfilt(d, c, v_acc);

x_acc = detrend(cumtrapz(t_acc, v_acc));
x_acc = filtfilt(d, c, x_acc);

end

