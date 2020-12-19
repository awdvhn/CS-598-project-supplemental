function [time, F_deNoise, x_acc, v_acc, a_acc, d, c,...
    iFstart, iFend, fc_LowPass] = computeDisp_AccRAWarr_R1_DetrendFilters(timeRAW, Fraw_Hammmer, aRAW_acc)

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
[f_FFTaccRAW, aFFT_accRAW] = regFFT(time, a_acc); 
[aFFT_Peaks, I_Peaks]  = findpeaks(aFFT_accRAW/max(aFFT_accRAW));
[aFFT_SortedPeaks, I_SortedPeaks] = sort(aFFT_Peaks);
aFFT_SortedPeaks = flip(aFFT_SortedPeaks);
I_SortedPeaks = flip(I_SortedPeaks);
fresMax_acc = max(f_FFTaccRAW(I_Peaks(I_SortedPeaks(aFFT_SortedPeaks >= 0.05))));
fresMin_acc = min(f_FFTaccRAW(I_Peaks(I_SortedPeaks(aFFT_SortedPeaks >= 0.5))));

% CusFigure('Frequency [Hz]',...
%     'Acceleration FFT [$\mathrm{m\cdot s^{-2}}$]', '');
% pCLR = distinguishable_colors(1);
% plot(f_FFTaccRAW, aFFT_accRAW, 'linewidth', 1.25, 'color', pCLR(1,:));
% xlim([0 1000])
% figure; plot(f_FFTacc(I_Peaks(I_SortedPeaks)), aFFT_SortedPeaks)

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
a_acc = filtfilt(b, a, a_acc); 

max1 = max(abs(detrend(aRAW_acc)));
if (max(abs(a_acc)) - max1)/max1 > 5e-3
    fprintf(['\t!!! Caution!!!\n\t Filtered Acceleration Max'...
        ' differs from  Raw Acceleration Max by more than: '...
        num2str(100*(max(abs(a_acc)) - max1)/max1) '%%\n'])
end

% %%% Plot of acceleration
% CusFigure('Time [s]',...
%     'Acceleration [$\mathrm{m\cdot s^{-2}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(time, detrend(aRAW_acc(I_ForceStart)), 'color', pCLR(1,:),...
%     'LineWidth', 1.25)
% plot(time, a_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
% legend('Raw', 'Filtered', 'box', 'off')
% xlim([0 1]); ylim(20*[-1 1])
% 
% xlim([0.15 0.5]); ylim(20*[-1 1])
% % 
% xlim([0 .1]); xticks(0:0.05:.1)
% ylim([-250 500]); yticks(-300:200:500)
% set(gcf, 'Position', [1 1 4 6])
% set(gca, 'Position', [.3 .2 .6 .7])
% 
%%% Plot of acceleration FFT
% [f_FFTacc, aFFT_acc] = regFFT(time, a_acc);
% CusFigure('Frequency [Hz]',...
%     'Acceleration FFT [$\mathrm{m\cdot s^{-2}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTaccRAW, aFFT_accRAW, 'linewidth', 1.25, 'color', pCLR(1,:));
% plot(f_FFTacc, aFFT_acc, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('Raw', 'Filtered', 'box', 'off')
% xlim([0 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Accelerometer Data -> Velocity 
fc_HighPass = 7.5;
Fc = fc_HighPass/Fnyq; % Filter cutoff ratio
[d,c] = butter(3,Fc,'high');
% %%% Filter Figure
% CusFigure('','','')
% freqz(d,c, 2e6, Fs); 
% xlim(fc_HighPass*[0 2]);
% ylim([40 180]); yticks(45:45:180)

v_accRAW = detrend(cumtrapz(time, a_acc));
% v_acc = detrend(cumtrapz(time, a_acc), 10);
v_acc = detrend(filtfilt(d, c, v_accRAW));

% %%% Plot of Velocity
% CusFigure('Time [s]',...
%     'Velocity [$\mathrm{m\cdot s^{-1}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(time, detrend(v_accRAW), 'color', pCLR(1,:),...
%     'LineWidth', 1.25)
% plot(time, v_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
% legend('Raw', 'Filtered', 'box', 'off')
% xlim([0 5]);
% 
% xlim([0 1]);
% 
% xlim([0 .01]);
% set(gcf, 'Position', [1 1 4 6])
% set(gca, 'Position', [.3 .2 .6 .7])
% 
% %%% Plot of Velocity FFT
% [f_FFTvel, aFFT_vel] = regFFT(time, v_acc);
% [f_FFTvelDet, aFFT_velDet] = regFFT(time, v_accDetrend);
% CusFigure('Frequency [Hz]',...
%     'Velocity FFT [$\mathrm{m\cdot s^{-1}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTvelDet, aFFT_velDet, 'linewidth', 1.25, 'color', pCLR(1,:));
% plot(f_FFTvel, aFFT_vel, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('Raw', 'Filtered', 'box', 'off')
% xlim([0 100])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Velocity -> Displacement 
fc_HighPass = 7.5;
Fc = fc_HighPass/Fnyq; % Filter cutoff ratio
[d,c] = butter(3,Fc,'high');
% %%% Filter Figure
% CusFigure('','','')
% freqz(d,c, 2e6, Fs); 
% xlim(fc_HighPass*[0 2]);
% ylim([40 180]); yticks(45:45:180)

x_accDetrend = detrend(cumtrapz(time, v_acc));
x_acc = filtfilt(d, c, x_accDetrend); 

%%% Plot of Displacement
CusFigure('Time [s]', 'Displacement [m]', '');
pCLR = distinguishable_colors(2);
plot(time, x_accDetrend, 'color', pCLR(1,:),...
    'LineWidth', 1.25)
plot(time, x_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
legend('Raw', 'Filtered', 'box', 'off')
xlim([0 5]);

xlim([0 1]); ylim(1.25e-3*[-1 1])
% 
% 
% %%% Plot of Displacement FFT
% [f_FFTdis, aFFT_dis] = regFFT(time, x_acc);
% [f_FFTdisDet, aFFT_disDet] = regFFT(time, x_accDetrend);
% CusFigure('Frequency [Hz]', 'Displacement FFT [m]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTdisDet, aFFT_disDet, 'linewidth', 1.25, 'color', pCLR(1,:));
% plot(f_FFTdis, aFFT_dis, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('Raw', 'Filtered', 'box', 'off')
% xlim([0 100])
% ylim([0 0.1e-3])
end

