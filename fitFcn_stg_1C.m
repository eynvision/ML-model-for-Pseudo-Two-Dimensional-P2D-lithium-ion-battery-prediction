function z = fitFcn_stg_ALL(x_sub,Simdata)
global Pstored fcnValStore count Crate

count = count+1;
fprintf('number of iterations: %d\n', count)

Pstored = [Pstored; x_sub] ;


load GME101_1C_EE_25oC_clean.mat;    
DC1=Data_exp;


Simdatac1 = MAIN_I_ROM_V3_1_1_PE(DC1,x_sub,1/1);



for i=1:length(Simdatac1.Vt)
    if isnan(Simdatac1.Vt(i))
       Simdatac1.Vt(i)=0;
    end
end

Fcnval = abs(real(sum((DC1.Voltage_V-Simdatac1.Vt).^2)));

fcnValStore = [fcnValStore; Fcnval];


z = real(Fcnval);
