function [y, State_est] = kalmanfilter(z, dt) 

% Initialize state transition matrix
A=[ 1         dt      0.5*dt^2 ;...     % [x  ]
    0         1         dt     ;...     % [Vx]
    0         0          1     ];     % [Ax]

H = [ 0 0 1 ];    % Initialize measurement matrix
Q = [0 0 0;...
     0 0 0;...
     0 0 1e-1];
R = 1e3 * eye(1);
persistent x_est p_est                % Initial state conditions
if isempty(x_est)
    x_est = zeros(3, 1);             % x_est=[x,y,Vx,Vy,Ax,Ay]'
    x_est(3) = z(1);
    p_est = zeros(3, 3);
end

% Predicted state and covariance
x_prd = A * x_est;
p_prd = A * p_est * A' + Q;

% Estimation
S = H * p_prd' * H' + R;
B = H * p_prd';
klm_gain = (S \ B)';

% Estimated state and covariance
x_est = x_prd + klm_gain * (z - H * x_prd);
p_est = p_prd - klm_gain * H * p_prd;

% Compute the estimated measurements
y = H * x_est;
State_est = x_est;
end                % of the function