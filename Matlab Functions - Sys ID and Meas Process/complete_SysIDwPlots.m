function [] = complete_SysIDwPlots(FileName)
%% Directories
inDir = ['1-Raw Measurements and First Processing/'...
    FileName '_Processed/MAT Files/']; % Directory in which experimental files are saved

outDir = ['3-System ID_CubicExponent_PositiveStiff/' FileName '/']; 
mkdir(outDir); % Directory where results are saved

%% Input Processed Data
time = load([inDir 'time.mat']); time = time.time;
F_ham = load([inDir 'F_ham.mat']); F_ham = F_ham.F_ham;
a = load([inDir 'a_acc_L2.mat']); a = a.a_acc_L2;
v = load([inDir 'v_acc_L2.mat']); v = v.v_acc_L2;
x = load([inDir 'x_acc_L2.mat']); x = x.x_acc_L2;

vmax = max(v(time > 0.1)); % Maximum Velocity

%%% Known Parameters
m = 2.097e-2; % [kg]
m_gdg = 18.85e-3; f_gdg = mean([repmat(24.17, [5,1]); repmat(24.19, [2,1]); repmat(24.25, [5,1])]);
k_g = m_gdg*(2*pi*f_gdg)^2; % [N/m]

zeta_gdg = mean([repmat(0.002704, [5,1]); repmat(0.003888, [2,1]); repmat(0.008665, [5,1])]);
d_g = m_gdg * 2 * zeta_gdg * (2*pi*f_gdg); % [Ns/m]

%% Stiffness Initial Guess

%%% Get Points for Restoring Effect: xdot = 0
fn_g = sqrt(k_g / m)/(2*pi);
t0 = time([find(time > 1/fn_g, 1, 'first'),...
    end - find(flip(x) > 1e-2*max(x), 1, 'first')]);
t0 = [0 1.75];
v0 = [0; 0]; % velocities are zero when considering the restoring effect only

Krest_0 = [k_g; 2.5e9; 3]; % parameters Bound
Krest_lb = [k_g; 1e7; 3];
Krest_ub = [10*k_g; 1e11; 3];

[ K_REST, R2_REST, P_REST, ti_REST, v_REST, x_REST, a_REST] = ...
    compute_InitialCoeffGuess( m, t0, v0, time, v, x, a,...
    @(Krest, m, x) aResFUN_Model(Krest, m, x), Krest_0, Krest_lb, Krest_ub);

%%% Nonlinear wire restoring effect fitting
a_RestWire = a_REST + k_g/m*x_REST;
a_RestWire = a_RestWire(5:end);

options = psoptimset;
options.Display = 'iter';
options.TolMesh = 1e-10;
options.MaxIter = 1e10;

Kwire0 = [0; 2.5e9; 3]; % parameters Bound
Kwire_lb = [0; 1e6; 3];
Kwire_ub = [10*k_g; 1e12; 3];

objFun = @(C) sum((a_RestWire - aResFUN_Model(C, m, x_REST(5:end))).^2)/sum((a_RestWire-mean(a_RestWire)).^2);

[K_Wire, Pres_acc] = patternsearch(objFun,...
    Kwire0,[],[],[],[],Kwire_lb,Kwire_ub,[],options);
R2_Wire = rsquare(a_RestWire, aResFUN_Model(K_Wire, m, x_REST(5:end)));

%% Damping Initial Guess 
% t0 = [0.7 2];
x0 = [0; 0]; % displacements are zero when considering the restoring effect only

aDampFUN_Model = @(C, m, v) -C(1)/m.*v; % functions

Krest_0 = d_g; % parameters Bound
Krest_lb = d_g;
Krest_ub = d_g*100;

[ Coefs_DAMP, R2_DAMP, P_DAMP, ti_DAMP, x_DAMP, v_DAMP , a_DAMP] = ...
    compute_InitialCoeffGuess( m, t0, x0, time, x, v, a,...
    aDampFUN_Model, Krest_0, Krest_lb, Krest_ub);

%% Plot Restoring Effect Points
close all;
pClrs = distinguishable_colors(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displacement Time Series
CusFigure('Time [s]', '$x$ [$\mathrm{m}$]','')
plot( time, x, '-', 'color', pClrs(1,:), 'linewidth', 1.5)
xlim(t0)

saveas(gcf, [outDir '4-x vs time Considered.fig'], 'fig')
saveas(gcf, [outDir '4-x vs time Considered.tiff'], 'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Points considered for Restoring Effect
tlim = [0.75 1];
plot_pointsForInitilGuess( time, a, v, x,...
    ti_REST, a_REST, v_REST, x_REST, tlim, vmax )
saveas(gcf, [outDir '5-Restoring Effect Points_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '5-Restoring Effect Points_Accelerometer.tiff'], 'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restoring Effect acceleration vs displacement
CusFigure('$x$ [$\mathrm{m}$]', '$a$ [$\mathrm{m\cdot s^{-2}}$]','Total Restoring Effect')

%%% Plot Fit
xResFit = linspace(1.1*min(x_REST), 1.1*max(x_REST), 1e3);
aResFit = aResFUN_Model(K_REST, m, xResFit);
p = plot(xResFit, aResFit, '-', 'color', pClrs(3,:), 'linewidth', 1.5);

%%% Plot Experimental
plot(x_REST, a_REST, 'x', 'color', pClrs(2,:), 'linewidth', 1.25, 'markersize', 10)

legend(p, {['$a=-\frac{\mathrm{' num2str(K_REST(1), '%.3E')...
    '}}{m}x-\frac{\mathrm{' num2str(K_REST(2), '%.3E')...
    '}}{m}x\cdot |x|^{\mathrm{' num2str(K_REST(3), '%i')...
    '-1}}$ ($\mathrm{R^2=' num2str(R2_REST, '%.3f') '}$)']},...
    'fontsize', 18,'location', 'northeast', 'Box', 'off')
saveas(gcf, [outDir '6a-Total Restoring Effect a vs x_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '6a-Total Restoring Effect a vs x_Accelerometer.tiff'], 'tiff')

%%% Restoring Effect of Nonlinear Wire
CusFigure('$x$ [$\mathrm{m}$]', '$a$ [$\mathrm{m\cdot s^{-2}}$]','Non-Linear Coupler Restoring Effect')
aFIT_Wire = aResFUN_Model(K_Wire, m, xResFit);
p = plot(xResFit, aFIT_Wire, '-', 'color', pClrs(3,:), 'linewidth', 1.5);
plot(x_REST(5:end), a_RestWire, 'x', 'color', pClrs(2,:), 'linewidth', 1.25, 'markersize', 10)

legend(p, {['$a=-\frac{\mathrm{' num2str(K_Wire(1), '%.3E')...
    '}}{m}x-\frac{\mathrm{' num2str(K_Wire(2), '%.3E')...
    '}}{m}x\cdot |x|^{\mathrm{' num2str(K_Wire(3), '%i')...
    '-1}}$ ($\mathrm{R^2=' num2str(R2_Wire, '%.3f') '}$)']},...
    'fontsize', 18,'location', 'northeast', 'Box', 'off')
saveas(gcf, [outDir '6b-NonLinear Coupler Restoring Effect a vs x_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '6b-NonLinear Coupler Restoring Effect a vs x_Accelerometer.tiff'], 'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Points considered for Damping Effect
tlim = [0.75 1];
plot_pointsForInitilGuess( time, a, v, x,...
    ti_DAMP, a_DAMP, v_DAMP, x_DAMP, tlim, vmax )
saveas(gcf, [outDir '7-Damping Effect Points_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '7-Damping Effect Points_Accelerometer.tiff'], 'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restoring Effect acceleration vs displacement
CusFigure('$v$ [$\mathrm{m\cdot s^{-1}}$]', '$a$ [$\mathrm{m\cdot s^{-2}}$]','')
plot( [-.15 .15], [0 0], 'color', [0.5 0.5 0.5])
plot( [0 0],[-2 1], 'color', [0.5 0.5 0.5])

%%% Plot Fit
vDampFit_acc = linspace(1.1*min(v_DAMP), 1.1*max(v_DAMP), 1e3);
aDampFit_acc = aDampFUN_Model(Coefs_DAMP, m, vDampFit_acc);
p = plot(vDampFit_acc, aDampFit_acc, '-', 'color', pClrs(3,:), 'linewidth', 1.5);

%%% Plot Data
plot(v_DAMP, a_DAMP, 'x', 'color', pClrs(2,:), 'linewidth', 1.25, 'markersize', 10)

legend(p, {['$a=-\frac{\mathrm{' num2str(Coefs_DAMP(1), '%.3E')...
    '}}{m}v$ ($\mathrm{R^2=' num2str(R2_DAMP, '%.3f') '}$)']},...
    'fontsize', 18,'location', 'southwest', 'Box', 'off')
saveas(gcf, [outDir '8a-Damping Effect a vs v_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '8a-Damping Effect a vs v_Accelerometer.tiff'], 'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot3(x, v, a)
grid on
xlabel('$x$ [m]')
ylabel('$v$ [$\mathrm{m\cdot s^{-1}}$]')
zlabel('$a$ [$\mathrm{m\cdot s^{-2}}$]')
saveas(gcf, [outDir '9-Restoring Force Curve_Accelerometer.fig'], 'fig')

%%
%% Nonlinear sytem-ID
close all

i_tend = find(time >= 0.5, 1, 'first');
i_Fend = find(abs(F_ham) > 0, 1, 'last');
t_Fend = time(i_Fend + 10);

n_downSample = 10;
time_fit = [time(1:i_Fend); downsample(time( (i_Fend+1):i_tend ), n_downSample)];
Fexp_fit = [F_ham(1:i_Fend); downsample(F_ham( (i_Fend+1):i_tend ), n_downSample)];
x_fit = [x(1:i_Fend); downsample(x( (i_Fend+1):i_tend ), n_downSample)];

% Objective Function
nonLinSys_Fun = @(t, X, Coefs, t_exp, F_exp, m) DynSys_SingleOsc_CubicNonlinear(t, X,...
    [Coefs(1), Coefs(2), Coefs(3), Coefs(4)], t_exp, F_exp, m, t_Fend);
objFun = @(C) objFun_2ndOrderSysID(C, time_fit, Fexp_fit, m, x_fit, nonLinSys_Fun);

%%% Parameters Bound
CnonLin_acc = [K_Wire(1) + k_g; K_Wire(2); K_Wire(3); d_g];% from zero crossings plot fitting
Krest_0 = CnonLin_acc;
Krest_lb = [k_g; CnonLin_acc(2)/5; 3; d_g];
Krest_ub = [2*CnonLin_acc(1); CnonLin_acc(2)*5; 3; 50*d_g];

%%% Global Optimizer
opt_Glob = optimoptions('surrogateopt', 'MaxFunctionEvaluations', 200);
prob_Glob = struct('objective',objFun, 'lb', Krest_lb, 'ub', Krest_ub,...
    'options', opt_Glob, 'solver','surrogateopt'); 
Cglob = surrogateopt(prob_Glob);
CnonLin_acc = Cglob';

%%% Local Optimizer
options = optimoptions('patternsearch');
options.Display = 'iter';
options.TolMesh = 1e-6;
options.MaxIter = 50e3;

for i_loc = 1:2
    Krest_0 = CnonLin_acc;
    Krest_lb = [k_g; CnonLin_acc(2)/5; 3; d_g];
    Krest_ub = [2*CnonLin_acc(1); CnonLin_acc(2)*5; 3; 50*d_g];
    [CnonLin_acc, PnonLin_acc] = patternsearch(objFun,...
        Krest_0,[],[],[],[],Krest_lb,Krest_ub,[],options);
end

coeff_Table = [Krest_lb, CnonLin_acc, Krest_ub]
save([outDir 'CoeffSummary_Accelerometer.mat'], 'coeff_Table')

%% Reconstruct the System
options2 = odeset; % ode45 Options
options2.RelTol = 1e-8;
options2.AbsTol = 1e-8;

X0 = [0; 0];
[~, Xrecon] = ode45(@(t,X) nonLinSys_Fun(t, X, CnonLin_acc,...
    time_fit, Fexp_fit, m), time_fit, X0, options2);
% [~, Xrecon] = ode45(@(t,X) DynSys_SingleOsc_CubicNonlinear(t, X, [CnonLin_acc(1), exp(CnonLin_acc(2)), CnonLin_acc(3), CnonLin_acc(4)],...
%     time_fit, Fexp_fit, m, t_Fend), time_fit, X0, options2);
xRecon = Xrecon(:,1);
R2_recon = rsquare(x_fit, xRecon);
T = table([Krest_lb; nan; nan], [CnonLin_acc; R2_recon; PnonLin_acc], [Krest_ub; nan; nan], 'VariableNames', {'Lower_Opt_Bd'; 'Opt_Result';'Upper_Opt_Bd'},...
    'RowNames',{'k_total [N/m]'; 'C [N/m3]'; 'beta'; 'd_total [Ns/m]'; 'R2'; 'Opt_ObjFun'})
writetable(T,[outDir 'Optimization_Results.csv'],'WriteRowNames',true)  

%%% Plot Reconstruction
CusFigure('Time [s]', '$x$ [m]',['$\mathrm{R^2='...
    num2str(R2_recon, '%.3f') '}$ - $v_{Max}=\mathrm{'...
    num2str(vmax, '%.1E') '\ m\cdot s^{-1}}$'])
plot(time, x, '-', 'color', pClrs(1,:), 'linewidth', 1.5)
plot(time_fit, xRecon, '-', 'color', pClrs(2,:), 'linewidth', 1.25)

xlim(time_fit([1 end]))
legend('Accelerometer', 'Reconstruction')

saveas(gcf, [outDir '10-Nonlinear Reconstruction_Accelerometer.fig'], 'fig')
saveas(gcf, [outDir '10-Nonlinear Reconstruction_Accelerometer.tiff'], 'tiff')

%%% Restoring Effect
CusFigure('$x$ [$\mathrm{m}$]', '$a$ [$\mathrm{m\cdot s^{-2}}$]','')

xResFit = linspace(1.1*min(x_REST), 1.1*max(x_REST), 1e3);
aResFit = aResFUN_Model(K_REST, m, xResFit);
p1 = plot(xResFit, aResFit, '-', 'color', pClrs(3,:), 'linewidth', 1.5);

xRes_ID = linspace(1.1*min(x_REST), 1.1*max(x_REST), 1e3)';
aRes_ID = aResFUN_Model([CnonLin_acc(1), CnonLin_acc(2), CnonLin_acc(3), CnonLin_acc(4)], m, xResFit);
p2 = plot(xRes_ID, aRes_ID, '-', 'color', pClrs(1,:), 'linewidth', 1.5);
for ii = 1:length(x_REST)
    [~, Ipts_4_r2IDrest(ii)] = min(abs(xRes_ID - x_REST(ii)));
end

plot(x_REST, a_REST, 'x', 'color', pClrs(2,:), 'linewidth', 1.25, 'markersize', 10)

legend([p1; p2], {['Direct Fit:' newline '$a=-\frac{\mathrm{' num2str(K_REST(1), '%.3E')...
    '}}{m}x-\frac{\mathrm{' num2str(K_REST(2), '%.3E')...
    '}}{m}x\cdot |x|^{\mathrm{' num2str(K_REST(3), '%i')...
    '-1}}$ ($\mathrm{R^2=' num2str(R2_REST, '%.3f') '}$)'];...
    ...
    [newline 'Reconstruction:' newline '$a=-\frac{\mathrm{' num2str(CnonLin_acc(1), '%.3E')...
    '}}{m}x-\frac{\mathrm{' num2str(CnonLin_acc(2), '%.3E')...
    '}}{m}x\cdot |x|^{\mathrm{' num2str(CnonLin_acc(3), '%i')...
    '-1}}$ ($\mathrm{R^2=' num2str(rsquare(a_REST, aRes_ID(Ipts_4_r2IDrest)'), '%.3f') '}$)']},...
    'fontsize', 16,'location', 'southwest', 'Box', 'off')
saveas(gcf, [outDir '11-Restoring Effect a vs x_with Reconst.fig'], 'fig')
saveas(gcf, [outDir '11-Restoring Effect a vs x_with Reconst.tiff'], 'tiff')

%%% Restoring Force Comparison for Identifed System
Frest_gd = - k_g*x;

k_WireLin = CnonLin_acc(1) - k_g;
Frest_wireLin = - k_WireLin*x;

C = CnonLin_acc(2);
beta = CnonLin_acc(3);
Frest_wireNonLin = - C*x.*abs(x).^(beta - 1);

CusFigure('x [m]', 'Restoring Force [N]', 'Identified System')
plot(x, Frest_gd, 'linewidth', 2, 'color', pClrs(1, :))
plot(x, Frest_wireLin, 'linewidth', 2, 'color', pClrs(2, :))
plot(x, Frest_wireNonLin, 'linewidth', 2, 'color', pClrs(3, :))
legend('Grounding Spring', 'Coupling Wire - Linear', 'Coupling Wire - Nonliear')
saveas(gcf, [outDir '12-Rest Force Dist F vs x_with Reconst.fig'], 'fig')
saveas(gcf, [outDir '12-Rest Force Dist F vs x_with Reconst.tiff'], 'tiff')

CusFigure('t [s]', 'Restoring [N]', 'Identified System')
plot(time, Frest_gd, 'linewidth', 2, 'color', pClrs(1, :))
plot(time, Frest_wireLin, '-', 'linewidth', 1.75, 'color', pClrs(2, :))
plot(time, Frest_wireNonLin, '-', 'linewidth', 1.5, 'color', pClrs(3, :))
legend('Grounding Spring', 'Coupling Wire - Linear', 'Coupling Wire - Nonliear')
saveas(gcf, [outDir '13-Rest Force Dist F vs time_with Reconst.fig'], 'fig')
saveas(gcf, [outDir '13-Rest Force Dist F vs time_with Reconst.tiff'], 'tiff')
end

