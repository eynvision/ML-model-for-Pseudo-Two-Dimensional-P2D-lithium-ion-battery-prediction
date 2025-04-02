function z = fitFcn_stg_ALL(x_sub,Simdata)
global Pstored fcnValStore count Crate

count = count+1;
fprintf('number of iterations: %d\n', count)

Pstored = [Pstored; x_sub] ;


load GME101_p3C_EE_25oC_discharging.mat; 
DC3=dchData;
load GME101_p2C_EE_25oC_discharging.mat;      
DC2=dchData;
load GME101_1C_EE_25oC_discharging.mat;    
DC1=dchData;
load GME101_1p5C_EE_25oC_discharging.mat;    
D2C=dchData;


Simdatac3 = MAIN_I_ROM_V3_1_1_PE(DC3,x_sub,1/3);
Simdatac2 = MAIN_I_ROM_V3_1_1_PE(DC2,x_sub,1/5);
Simdatac1 = MAIN_I_ROM_V3_1_1_PE(DC1,x_sub,1/1);
Simdata2c = MAIN_I_ROM_V3_1_1_PE(D2C,x_sub,1.5/1);

for i=1:length(Simdatac1.Vt)
    if isnan(Simdatac1.Vt(i))
       Simdatac1.Vt(i)=0;
    end
end
for i=1:length(Simdatac2.Vt)
    if isnan(Simdatac2.Vt(i))
       Simdatac2.Vt(i)=0;
    end
end
for i=1:length(Simdatac3.Vt)
    if isnan(Simdatac3.Vt(i))
       Simdatac3.Vt(i)=0;
    end
end
for i=1:length(Simdata2c.Vt)
    if isnan(Simdata2c.Vt(i))
       Simdata2c.Vt(i)=0;
    end
end
Fcnval = abs(real(sum((DC3.Voltage-Simdatac3.Vt).^2)))+abs(real(sum((DC2.Voltage-Simdatac2.Vt).^2)))+abs(real(sum((DC1.Voltage-Simdatac1.Vt).^2)))+abs(real(sum((D2C.Voltage-Simdata2c.Vt).^2)));

fcnValStore = [fcnValStore; Fcnval];


z = real(Fcnval);
