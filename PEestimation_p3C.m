clear; clc; close all;
tic
warning on
%% ParaID
global Pstored  fcnValStore count  xglobal Crate

n1=10; %number of veriable

Pstored = zeros(1,n1);
fcnValStore = [];
count = 0;

opts = gaoptimset('generations',35,'TolFun',0.05);

range1=[0.1e-10    0.1e-10   0.1e-6    0.40  1e-5  1e-5 90  0.45 100e-6 100e-6];

range2=[10e-9    10e-9     10e-4      0.65   90e-4 50e-4 118.62  0.65 1000e-5 1000e-5];


Crate_list=[1/3];

tic
for Crate=Crate_list
    [xglobal] = ga(@fitFcn_stg_p3C, n1, [],[],[],[],range1,range2,[],[],opts);
    switch Crate
        case 1/3
            save(['PM_EE_p3C_new.mat'], 'xglobal','Pstored', 'fcnValStore', 'count');
            load('GME101_p3C_EE_25oC_clean.mat')
            Simdata = MAIN_I_ROM_V3_1_1_PE(Data_exp,xglobal,Crate);
            figure; plot(Data_exp.time/60,Data_exp.Vt); hold on; plot(Simdata.t/60,Simdata.Vt);
            plot(Simdata.t/60,Simdata.OCV); legend('V_{exp}','V_{sim}','OCV'); title('C/3');
            figure; plot(Data_exp.time/60,Data_exp.Vt-Simdata.Vt);
        
    end
end
toc

