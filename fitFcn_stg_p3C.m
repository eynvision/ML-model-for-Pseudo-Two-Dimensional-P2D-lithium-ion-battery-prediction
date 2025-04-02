function z = fitFcn_stg_ALL(x_sub,Simdata)
global Pstored fcnValStore count Crate

count = count+1;
fprintf('number of iterations: %d\n', count)

Pstored = [Pstored; x_sub] ;


load GME101_p3C_EE_25oC_clean.mat; 
DC3=Data_exp;



Simdatac3 = MAIN_I_ROM_V3_1_1_PE(DC3,x_sub,1/3);


for i=1:length(Simdatac3.Vt)
    if isnan(Simdatac3.Vt(i))
       Simdatac3.Vt(i)=0;
    end
end

Fcnval = abs(real(sum((DC3.Voltage_V-Simdatac3.Vt).^2)));

fcnValStore = [fcnValStore; Fcnval];


z = real(Fcnval);
