
clear all, close all; clc;

%% ===================== 1) C/3 section =====================
Crate = 1/5;   % (User-labeled as '1/3C' in your code)
load('PM_EE_p2C_new.mat');
% xglobal(2) = 5.27782380580800e-09;
%9.049760877112130e-09
% xglobal(6) = xglobal(6) - 0.00009;

load('GME101_p2C_EE_25oC_clean.mat');
Simdata = MAIN_I_ROM_V3_1_1_PE(Data_exp, xglobal, Crate);

% Compute energies for bar chart
W_in_C3      = trapz(Simdata.t, Simdata.Vt .* abs(Simdata.I));
W_heat_irr_C3= trapz(Simdata.t, Simdata.Heat_irr_tot);

figure; 
plot(Data_exp.Time_S(1:end-10), Data_exp.Voltage_V(1:end-10),'r--'); 
hold on; 
plot(Simdata.t(1:end-10), Simdata.Vt(1:end-10),'r');
Plot_default;

figure(2); 
plot(Data_exp.Time_S(1:end-10), Data_exp.Total_HGR_W(1:end-10),'r--'); 
hold on; 
plot(Simdata.t(1:end-10), Simdata.Heat_tot(1:end-10),'r');
Plot_default;


%% ===================== 2) C/2 section =====================
Crate = 1/3;   % (User-labeled as '1/2C' in your code)
load('PM_EE_p3C_new.mat');

load('GME101_p3C_EE_25oC_clean.mat')
Simdata = MAIN_I_ROM_V3_1_1_PE(Data_exp, xglobal, Crate);

% Compute energies
W_in_C2       = trapz(Simdata.t, Simdata.Vt .* abs(Simdata.I));
W_heat_irr_C2 = trapz(Simdata.t, Simdata.Heat_irr_tot);

figure(1);
plot(Data_exp.Time_S(1:end-10), Data_exp.Voltage_V(1:end-10),'g--'); 
hold on;
plot(Simdata.t(1:end-10), Simdata.Vt(1:end-10),'g');
Plot_default;

figure(2);
plot(Data_exp.Time_S(1:end-10), Data_exp.Total_HGR_W(1:end-10),'g--'); 
hold on;
plot(Simdata.t(1:end-10), Simdata.Heat_tot(1:end-10),'g');
Plot_default;


%% ===================== 3) 1C section =====================
Crate = 1;
load('PM_EE_ALL_new.mat');

load('GME101_1C_EE_25oC_clean.mat')

Simdata = MAIN_I_ROM_V3_1_1_PE(Data_exp, xglobal, Crate);

% Compute energies
W_in_1C       = trapz(Simdata.t, Simdata.Vt .* abs(Simdata.I));
W_heat_irr_1C = trapz(Simdata.t, Simdata.Heat_irr_tot);

figure(1);
plot(Data_exp.Time_S(1:end-10), Data_exp.Voltage_V(1:end-10),'m--'); 
hold on; 
plot(Simdata.t(1:end-10), Simdata.Vt(1:end-10),'m');
Plot_default;

figure(2);
plot(Data_exp.Time_S(1:end-10), Data_exp.Total_HGR_W(1:end-10),'m--'); 
hold on; 
plot(Simdata.t(1:end-10), Simdata.Heat_tot(1:end-10),'m');
Plot_default;
