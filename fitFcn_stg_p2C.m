function z = fitFcn_stg_ALL(x_sub,Simdata)
global Pstored fcnValStore count Crate

count = count+1;
fprintf('number of iterations: %d\n', count)

Pstored = [Pstored; x_sub] ;



load GME101_p2C_EE_25oC_clean.mat;      
DC2=Data_exp;


Simdatac2 = MAIN_I_ROM_V3_1_1_PE(DC2,x_sub,1/5);


for i=1:length(Simdatac2.Vt)
    if isnan(Simdatac2.Vt(i))
       Simdatac2.Vt(i)=0;
    end
end

Fcnval = abs(real(sum((DC2.Voltage_V-Simdatac2.Vt).^2)));

fcnValStore = [fcnValStore; Fcnval];


z = real(Fcnval);
