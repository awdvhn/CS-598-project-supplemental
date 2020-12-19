function errorSUM = objFun_2ndOrderSysID(C, t_exp, F_exp, m, x_exp, sysFUN)
% Run the system simulation
options2 = odeset; % ode45 Options
options2.RelTol = 1e-8;
options2.AbsTol = 1e-8;

X0 = [0; 0];
[~, Xsim] = ode45(@(t,X) sysFUN(t, X, C, t_exp, F_exp, m), t_exp, X0, options2);

x_sim = Xsim(:,1);

errorSUM = sum((x_exp - x_sim).^2)./sum((x_exp - mean(x_exp)).^2);

end