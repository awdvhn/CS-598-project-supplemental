function [time, x_vib, v_vib] = computeDisp_VibRAWarr(timeRAW, vRAW_vib, iFstart)

%% Trimming to Range of interest
%%% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

%%% Trim Time Instants to after Hit Region
time = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

%%% Trim Accelerometer data to after Hit Region
v_vib = detrend(vRAW_vib(I_ForceStart));

%% Create Filter for Processing
Fs = 1/mean(diff(time));
Fnyq = Fs/2;

%%% Approximate of natural frequency
[f_FFTvib, vFFT_vib] = regFFT(time, v_vib); 
[vFFT_Peaks, I_Peaks]  = findpeaks(vFFT_vib/max(vFFT_vib));
[vFFT_SortedPeaks, I_SortedPeaks] = sort(vFFT_Peaks);
vFFT_SortedPeaks = flip(vFFT_SortedPeaks);
I_SortedPeaks = flip(I_SortedPeaks);
fresMax_acc = max(f_FFTvib(I_Peaks(I_SortedPeaks(vFFT_SortedPeaks >= 0.05))));
fresMin_acc = min(f_FFTvib(I_Peaks(I_SortedPeaks(vFFT_SortedPeaks >= 0.5))));

%%% Low-Pass Filter for Experimental Noise
fc_LowPass = 12500;
if fc_LowPass < 5*fresMax_acc
    fprintf(['\t !!! Caution !!!\n\t Low Pass Filter: {Fc = ' num2str(fc_LowPass)...
        ' Hz} < {5*Important Dynamics Frequency = ' num2str(5*fresMax_acc) 'Hz}\n'])
end
Fc_LowPass = fc_LowPass/Fnyq;
[b,a] = butter(11,Fc_LowPass,'low');
% figure; freqz(b, a, 1e4, Fs); xlim([0 2*fc_LowPass]); ylim([-10 5])

%%% High-Pass Filter for Integration
fc_HighPass = 5;
if fc_HighPass > fresMin_acc/3
    fprintf(['\t !!! Caution !!!\n\t High Pass Filter: {Fc = ' num2str(fc_HighPass)...
        ' Hz} > {5*Important Dynamics Frequency = ' num2str(fresMin_acc/3) 'Hz}\n'])
end
Fc = fc_HighPass/Fnyq; % Filter cutoff ratio
[d,c] = butter(2,Fc,'high');
% freqz(d,c, 2e6, Fs); xlim(fc_HighPass*[0 2]); ylim([-10 5])

%% Integrate to get displacement

%%% Velocity Data -> Displacement 
v_vib = filtfilt(b, a, v_vib); 
max1 = max(abs(detrend(vRAW_vib)));
if (max(abs(v_vib)) - max1)/max1 > 5e-3
    fprintf(['\t!!! Caution!!!\n\t Filtered Velocity Maximum differs from '...
        'Raw Velocity Maximum by more than: ' num2str(100*(max(abs(v_vib)) - max1)/max1) '%%\n'])
end
% figure;plot(timeRAW -timeRAW(iFstart), detrend(vRAW_vib), time, v_vib)

x_vib = detrend(cumtrapz(time, v_vib));
x_vib = filtfilt(d, c, x_vib);
% figure; plot(time, x_vib); grid on

end

