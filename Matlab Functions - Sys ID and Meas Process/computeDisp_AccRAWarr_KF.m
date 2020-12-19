function [time, F_deNoise, x_acc, v_acc, a_acc,...
    iFstart, iFend] = computeDisp_AccRAWarr_KF(timeRAW, Fraw_Hammmer, aRAW_acc)

%%% Denoise Hammer Force, Get Impulse Start and End instants after denoise
[F_deNoise, iFstart, iFend] = hammerForce_trimDeNoise(Fraw_Hammmer); 

%%% Index Range after the hit 
I_ForceStart = iFstart:length(timeRAW);

%%% Trim Time Instants to after Hit Region
time = timeRAW(I_ForceStart) - timeRAW(I_ForceStart(1));

%%% Trim Denoised Force to after Hit Region
F_deNoise = F_deNoise(I_ForceStart);

%%% Trim Accelerometer data to after Hit Region
a_acc_Detrend = detrend(aRAW_acc(I_ForceStart));

%% Apply Kalman Filter
%  Remove acceleration noiseand estimate velocity & displacement
dt = mean(diff(time));
Z = a_acc_Detrend';
S_est = zeros(3, length(a_acc_Detrend));
for ii_t = 1:length(a_acc_Detrend)
    [~, S_ii] = kalmanfilter(Z(:,ii_t), dt);
    S_est(:, ii_t) = S_ii;
end
a_acc = S_est(3,:);
v_acc = S_est(2,:);
x_acc = S_est(1,:);

CusFigure('Time [s]', 'Acceleration [m/s2]', '')
plot(time, a_acc_Detrend)
plot(time, a_acc)
xlim([0 1])
ylim([-50 50])

CusFigure('Time [s]', 'Velocity [m/s]', '')
plot(time, v_acc)
% xlim([0 1])
% ylim([-50 50])

CusFigure('Time [s]', 'Displacement [m/s]', '')
plot(time, detrend(x_acc, 0))
end

