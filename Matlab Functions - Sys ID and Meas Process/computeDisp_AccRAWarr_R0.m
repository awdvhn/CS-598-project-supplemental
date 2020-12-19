function [time, F_deNoise, x_acc, v_acc, a_acc, d, c,...
    iFstart, iFend, fc_LowPass] = computeDisp_AccRAWarr_R0(timeRAW, Fraw_Hammmer, aRAW_acc)

%%% Denoise Hammer Force, Get Impulse Start and End instants after denoise
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(Fraw_Hammmer); 

%%% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

%%% Trim Time Instants to after Hit Region
time = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

%%% Trim Denoised Force to after Hit Region
F_deNoise = F_deNoise(I_ForceStart);

%%% Trim Accelerometer data to after Hit Region
a_acc = detrend(aRAW_acc(I_ForceStart));

%% Create Filter for Processing
Fs = 1/mean(diff(time));
Fnyq = Fs/2;

%%% Approximate of natural frequency
[f_FFTacc, aFFT_acc] = regFFT(time, a_acc); 
[aFFT_Peaks, I_Peaks]  = findpeaks(aFFT_acc/max(aFFT_acc));
[aFFT_SortedPeaks, I_SortedPeaks] = sort(aFFT_Peaks);
aFFT_SortedPeaks = flip(aFFT_SortedPeaks);
I_SortedPeaks = flip(I_SortedPeaks);
fresMax_acc = max(f_FFTacc(I_Peaks(I_SortedPeaks(aFFT_SortedPeaks >= 0.05))));
fresMin_acc = min(f_FFTacc(I_Peaks(I_SortedPeaks(aFFT_SortedPeaks >= 0.5))));

% figure; plot(f_FFTacc, aFFT_acc)
% figure; plot(f_FFTacc(I_Peaks(I_SortedPeaks)), aFFT_SortedPeaks)

%%% Band pass filter to eliminate accelerometer related frequency
f1 = 350/Fnyq;
f2 = 2300/Fnyq;
[f,e] = butter(3,[f1 f2],'stop');
% figure; freqz(f, e, 1e4, Fs); xlim([0 2*f2]*Fnyq); ylim([-10 5])

%%% Low-Pass Filter for Experimental Noise
% fc_LowPass = 12500;
fc_LowPass = 250;
if fc_LowPass < 5*fresMax_acc
    fprintf(['\t !!! Caution !!!\n\t Low Pass Filter: {Fc = ' num2str(fc_LowPass)...
        ' Hz} < {5*Important Dynamics Frequency = ' num2str(5*fresMax_acc) 'Hz}\n'])
end
Fc_LowPass = fc_LowPass/Fnyq;
[b,a] = butter(7,Fc_LowPass,'low');
% figure; freqz(b, a, 1e4, Fs); xlim([0 2*fc_LowPass]); ylim([-10 5])

%% Integrate to get displacement

%%% Accelerometer Data -> Displacement 
% a_acc = filtfilt(b, a, a_acc);

[~, i_MidForce] = max(F_deNoise);
a_toBeFiltered = a_acc(i_MidForce:end);
i_a0_First = find(sign(a_toBeFiltered) ~= sign(a_toBeFiltered(1)), 1, 'first');
a_toBeFiltered = a_toBeFiltered(i_a0_First:end);

a_toBeFiltered = filtfilt(b, a, a_toBeFiltered); 
a_Forced = a_acc(1:(i_MidForce+i_a0_First-2)); 
t_Forced = time(1:(i_MidForce+i_a0_First-2));
a_Forced = smooth(t_Forced, a_Forced, floor((iFend-iFstart)/2),'loess');
a_Forced(a_Forced < a_toBeFiltered(1)) = a_toBeFiltered(1);
a_acc = [a_Forced; a_toBeFiltered];

% a_Forced = a_acc(1:(iFend-iFstart+1)); 
% t_Forced = time(1:(iFend-iFstart+1));
% a_Forced = smooth(t_Forced, a_Forced, floor((iFend-iFstart)/3),'loess');
% 
% a_notForced = a_acc((iFend-iFstart+2):end);
% t_notForced = time((iFend-iFstart+2):end);
% a_notForced = smooth(t_notForced, a_notForced, ceil((iFend-iFstart)*20),'loess');
% a_Forced(a_Forced < 0 ) = a_notForced(1);
% 
% a_acc = [a_Forced; a_notForced];
% % a_acc = movmean(a_acc, ceil(*(iFend-iFstart)));

%%% Predict inital Acceleration from F = ma
a_toBeFiltered = a_acc((iFend-iFstart+2):end);
a_toBeFiltered = filtfilt(b, a, a_toBeFiltered); 
a_Forced = F_deNoise(1:(iFend-iFstart+1))/.0245; 
a_Forced(a_Forced < a_toBeFiltered(1)) = a_toBeFiltered(1);
a_acc = [a_Forced; a_toBeFiltered];

max1 = max(abs(detrend(aRAW_acc)));
if (max(abs(a_acc)) - max1)/max1 > 5e-3
    fprintf(['\t!!! Caution!!!\n\t Filtered Acceleration Max differs from '...
        'Raw Acceleration Max by more than: ' num2str(100*(max(abs(a_acc)) - max1)/max1) '%%\n'])
end
% figure;plot(timeRAW -timeRAW(iFstart), detrend(aRAW_acc), time, a_acc); grid on; xlim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Accelerometer Data -> Velocity 
fc_HighPass = 13;
Fc = fc_HighPass/Fnyq; % Filter cutoff ratio
[d,c] = butter(4,Fc,'high');
% freqz(d,c, 2e6, Fs); xlim(fc_HighPass*[0 2]); ylim([-10 5])
v_acc = detrend(cumtrapz(time, a_acc), 5);
% v_acc = detrend(filtfilt(d, c, v_acc));

% figure; plot(time, detrend(cumtrapz(time, a_acc)), time, v_acc); grid on; xlim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Velocity -> Displacement 
fc_HighPass = 3;
Fc = fc_HighPass/Fnyq; % Filter cutoff ratio
[d,c] = butter(4,Fc,'high');
% freqz(d,c, 2e6, Fs); xlim(12.5*[0 2]); ylim([-10 5])
x_acc = detrend(cumtrapz(time, v_acc));
x_acc = filtfilt(d, c, x_acc);

% figure; plot(time,detrend(cumtrapz(time, v_acc)), time, x_acc); grid on; xlim([0 1])

end

