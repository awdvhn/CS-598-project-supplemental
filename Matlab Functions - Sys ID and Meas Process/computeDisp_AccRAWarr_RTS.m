function [time, F_deNoise, x_acc, v_acc, a_acc, iFstart, iFend]...
    = computeDisp_AccRAWarr_RTS(timeRAW, Fraw_Hammmer, aRAW_acc)

%%% Denoise Hammer Force, Get Impulse Start and End instants after denoise
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(Fraw_Hammmer); 

%%% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

%%% Trim Time Instants to after Hit Region
time = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

%%% Trim Denoised Force to after Hit Region
F_deNoise = F_deNoise(I_ForceStart);

%% Low-Pass Filter to compute Low Frequency Components in RTS
a_acc_Detrend = detrend(aRAW_acc(I_ForceStart));
Fs = 1/mean(diff(time));
Fnyq = Fs/2;

fc_veryLowPass = 2.5;
Fc_veryLowPass = fc_veryLowPass/Fnyq;
[A,B,C,D] = butter(1,Fc_veryLowPass,'low');

% %%% Generate Figure to Double Check State-Space Filter Applied
% [b, a] = butter(1,Fc_veryLowPass,'low');
% X_Filter = zeros(size(a_acc_Detrend) + 1);
% a_Filtered = nan*ones(size(a_acc_Detrend));
% for ii_acc = 1:length(a_acc_Detrend)
%     X_Filter(ii_acc + 1) = A*X_Filter(ii_acc) + B*a_acc_Detrend(ii_acc);
%     a_Filtered(ii_acc) = C*X_Filter(ii_acc) + D*a_acc_Detrend(ii_acc);
% end
% CusFigure('Time [s]', 'Acceleration [$\mathrm{m\cdot s^{-2}}$]', '')
% pCLR = distinguishable_colors(3);
% plot(time, a_acc_Detrend, 'linewidth', 1.25, 'color', pCLR(1, :) )
% plot(time, a_Filtered, 'linewidth', 1.25, 'color', pCLR(2, :) )
% plot(time, filter(b, a, a_acc_Detrend), ':', 'linewidth', 2.25,...
%     'color', pCLR(3, :) )
% legend('Raw', 'State-Space Low Freq. Comp.',...
%     'Transfer Function Low Freq. Comp.', 'box', 'off')
% ylim(20*[-1 1])
% xlim([0 1])

%% Apply Rauch-Tung-Striebel Smoothing 
addpath(genpath('KalmanAll'))

%%% Define Smoother System 
% X(t+1) = F X(t) + noise(Q) 
% X(t) = [x(t); v(t); a(t); x_lowFreq(t); v_lowFreq(t); a_lowFreq(t)]
% Y(t) = H X(t) + noise(R) with Y(t) = Measured acceleration
y = [a_acc_Detrend'; zeros(2, length(a_acc_Detrend))];
dt = mean(diff(time));
F = [ 1,     dt,   0.5*dt^2, 0, 0;... % [x] % The system matrix
      0,     1,       dt,    0, 0;... % [v]
      0,     0,       1,     0, 0;... % [a]
      B,     0,       0,     A, 0;... % [x_lowFreq(k+1)]
      0,     B,       0,     0, A] ;... % [v_lowFreq(k+1)]
%       0,     0,       B,     0, 0, A];    % [a_lowFreq(k+1)]

H = [0, 0, 1, 0, 0;... % Measured acceleration(k)
     D, 0, 0, C, 0;... % x_lowFreq(k)
     0, D, 0, 0, C];... % v_lowFreq(k)
%      0, 0, D, 0, 0, C];    % a_lowFreq(k)
     
ss = size(F, 2); % State space size
os = size(y, 1); % observation space size
Q = 2e-4*eye(ss); % The system covariance
R = 3e-1*eye(os); % The observation covariance
X0 = zeros(ss, 1); % the initial state (column) vector
cov0 = 8.2e2*eye(ss); % the initial state covariance

%%% Apply Smoother
S_smooth = kalman_smoother(y, F, H, Q, R, X0, cov0);
rmpath(genpath('KalmanAll'))

%%% Assign to output variables
x_acc = S_smooth(1,:);
v_acc = S_smooth(2,:);
a_acc = S_smooth(3,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot of Acceleration
% CusFigure('Time [s]',...
%     'Acceleration [$\mathrm{m\cdot s^{-2}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(time, a_acc_Detrend, 'linewidth', 1.25, 'color', pCLR(1,:),...
%     'LineWidth', 1.25)
% plot(time, a_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
% legend('Raw', 'RTS', 'box', 'off')
% xlim([0 1]); ylim(20*[-1 1])
% 
% xlim([0.15 0.5]); ylim(20*[-1 1])
% 
% xlim([0 .1]);xticks(0:0.05:0.1);
% ylim([-250 500]); yticks(-200:200:500)
% set(gcf, 'Position', [1 1 4 6])
% set(gca, 'Position', [.3 .2 .6 .7])
% 
% %%% Plot of acceleration FFT
% [f_FFTacc, aFFT_acc] = regFFT(time, a_acc);
% [f_FFTaccRAW, aFFT_accRAW] = regFFT(time, a_acc_Detrend);
% CusFigure('Frequency [Hz]',...
%     'Acceleration FFT [$\mathrm{m\cdot s^{-2}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTaccRAW, aFFT_accRAW, 'linewidth', 1.25, 'color', pCLR(1,:));
% plot(f_FFTacc, aFFT_acc, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('Raw', 'RTS', 'box', 'off')
% xlim([0 1000])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot of Velocity
% CusFigure('Time [s]',...
%     'Velocity [$\mathrm{m\cdot s^{-1}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(time, v_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
% legend('RTS', 'box', 'off')
% xlim([0 5]);
% 
% xlim([0 1]);
% 
% xlim([0 .1]);xticks(0:0.05:0.1);
% set(gcf, 'Position', [1 1 4 6])
% set(gca, 'Position', [.3 .2 .6 .7])
% 
% %%% Plot of Velocity FFT
% [f_FFTvel, aFFT_vel] = regFFT(time, v_acc);
% CusFigure('Frequency [Hz]',...
%     'Velocity FFT [$\mathrm{m\cdot s^{-1}}$]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTvel, aFFT_vel, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('RTS', 'box', 'off')
% xlim([0 100])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot of Displacement
% CusFigure('Time [s]', 'Displacement [m]', '');
% pCLR = distinguishable_colors(2);
% plot(time, x_acc, 'color', pCLR(2,:), 'LineWidth', 1.25)
% legend('RTS', 'box', 'off')
% xlim([0 5]);
% 
% xlim([0 1]); ylim(1.25e-3*[-1 1])
% 
% 
% %%% Plot of Displacement FFT
% [f_FFTdis, aFFT_dis] = regFFT(time, x_acc);
% CusFigure('Frequency [Hz]', 'Displacement FFT [m]', '');
% pCLR = distinguishable_colors(2);
% plot(f_FFTdis, aFFT_dis, 'linewidth', 1.25, 'color', pCLR(2,:));
% legend('RTS', 'box', 'off')
% xlim([0 100])
% ylim([0 0.1e-3])
end

