clear
close all
format long

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');

addpath(genpath('Matlab Functions - Figures'));
addpath(genpath('Matlab Functions - Sys ID and Meas Process'));
clc

%%
inParentDir = '1-Matlab Exp and RTS Results/';
outdir = '2-Sim for Exp Sim/';
mkdir(outdir)

%%% Load Forcing Measurement
load([inParentDir 'F_ham.mat'])
load([inParentDir 'time.mat'])

%%% Simulation Time [s]
T_sim = 1;
N = 7;
loc = [1; 0; 0; 0];

%% Assign system with Experimental Parameters

%%% Masses
m1 = 0.0245 * ones([7 1]);
m2 = 0.0245 * ones([7 1]);

%%% Grounding Springs
kg1 = [472.61 1791.9 1878.8 1800.1 1772.6 1841.5 1794.2]';
kg2 = [487.29 495.92 494.58 429.6 419.5 495.92 466.6]';
d_g1 = [0.0177 0.0550 0.0556 0.0255 0.0255 0.0213 0.0326]';
d_g2 = [0.0134 0.0322 0.0235 0.0214 0.0182 0.0322 0.0165]';

%%% Linear Couplers
d_c = [0 0 0 0 0 0 0]';
coupling=[115.4 131.9 110.8 112.3 124.8 124.4 128.6]';

%%% Nonlinear Coupler
beta = 2.902;
d_NL1 = 1.61/2 * 0.02 * ones([6 1]);
d_NL2 = 1.61/2 * 0.02 * ones([6 1]);
kNL1 = 2.3e8 * ones([6 1]);
kNL2 = 2.3e8 * ones([6 1]);
kNL_L1 = 0 * -24.4 * ones([6 1]);
kNL_L2 = 0 * -24.4 * ones([6 1]);

%% Simulate the Experimental System
i_tend = find(time >= 0.01, 1, 'first');
i_send = find(time >= T_sim, 1, 'first');
Tspan1=time(1:i_tend);
Tspan2=time(i_tend:i_send);

y0=zeros(4*N,1);
v0=trapz(Tspan1,F_ham(1:i_tend))/m1(1);
options=odeset('Reltol',1e-10,'Abstol',1e-10);
[T1,Y1]=ode45(@(t,y) dydtarbi(t,y,N,kg1,kg2,d_g1,d_g2,d_c,d_NL1,d_NL2,m1,m2,coupling,kNL1,kNL2,kNL_L1,kNL_L2,beta,F_ham(1:i_tend),Tspan1, loc),Tspan1,y0,options);
y_intermediate=Y1(end,:)';
[T2,Y2]=ode45(@(t,y) dydtunforced(t,y,N,kg1,kg2,d_g1,d_g2,d_c,d_NL1,d_NL2,m1,m2,coupling,kNL1,kNL2,kNL_L1,kNL_L2,beta),Tspan2,y_intermediate,options);
T=[T1(1:end-1);T2];
Y=[Y1(1:end-1,:);Y2];

%% Plot Figure
pClrs = distinguishable_colors(14);

figure
h17=plot(T,Y(:,25)-12*v0/150,'color',pClrs(14,:));
hold on
h16=plot(T,Y(:,21)-10*v0/150,'color',pClrs(13,:));
hold on
h15=plot(T,Y(:,17)-8*v0/150,'color',pClrs(12,:));
hold on
h14=plot(T,Y(:,13)-6*v0/150,'color',pClrs(11,:));
hold on
h13=plot(T,Y(:,9)-4*v0/150,'color',pClrs(10,:));
hold on
h12=plot(T,Y(:,5)-2*v0/150,'color',pClrs(9,:));
hold on
h11=plot(T,Y(:,1),'color',pClrs(8,:));
hold on
h21=plot(T,Y(:,3)+2*v0/150,'color',pClrs(7,:));
hold on
h22=plot(T,Y(:,7)+4*v0/150,'color',pClrs(6,:));
hold on
h23=plot(T,Y(:,11)+6*v0/150,'color',pClrs(5,:));
hold on
h24=plot(T,Y(:,15)+8*v0/150,'color',pClrs(4,:));
hold on
h25=plot(T,Y(:,19)+10*v0/150,'color',pClrs(3,:));
hold on
h26=plot(T,Y(:,23)+12*v0/150,'color',pClrs(2,:));
hold on
h27=plot(T,Y(:,27)+14*v0/150,'color',pClrs(1,:));
hold off
set(gca,'fontsize',16);
axis([0 max(T) -13*v0/150 15*v0/150]);
%str=['\zeta=',num2str(xi/2)];
title(['Simulation'])
xlim([0 1])
xlabel('Time (s)','Fontsize',20);
ylabel('Displacement (m)','Fontsize',20);
legend([h27 h26 h25 h24 h23 h22 h21 h11 h12 h13 h14 h15 h16 h17],'L7','L6','L5','L4','L3','L2','L1','R0','R1','R2','R3','R4','R5','R6')
set(gca, 'Units', 'inches', 'Position', [1.25 0.75 2.5 6])
set(gcf,'Units', 'inches', 'Position', [0.75 0.75 5 7.25])
legend('units', 'inches', 'position',  [4 0.725 0.5 6.05])
saveas(gcf,[outdir,'time_series_simulation.fig'])
% close all

%% Save Data for Regression

%%% Extract Velocity and Displacement of R0
xSim_R0 = Y(:,1);
vSim_R0 = Y(:,2);

%%% Compute Accelertion 
aSim_R0 = loc(1).*F_ham(1:i_send)/m1(1) + (-kg1(1).*Y(:, 1)-kNL1(1).*...
    abs(Y(:, 1)-Y(:, 5)).^(beta-1).*(Y(:, 1)-Y(:, 5))-kNL_L1(1).*(Y(:, 1)-Y(:, 5))...
    -coupling(1)*(Y(:, 1)-Y(:, 3))-d_g1(1)*Y(:, 2)-d_c(1)*(Y(:, 2)-Y(:, 4))...
    -d_NL1(1).*(Y(:, 2)-Y(:, 6)))/m1(1);

%%% Save dat
TFAVX = [T, F_ham(1:i_send), aSim_R0, vSim_R0, xSim_R0];
writematrix(TFAVX,[outdir, 'SimTFAVX_R0.txt']);
save([outdir, 'SimMatTFAVX_R0'], 'TFAVX');

%%
rmpath(genpath('Functions - Figures'));
rmpath(genpath('Functions - Sys ID and Meas Process'));

