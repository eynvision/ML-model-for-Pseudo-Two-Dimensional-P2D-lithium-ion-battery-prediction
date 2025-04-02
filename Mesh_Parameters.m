function [PM]= Mesh_Parameters( PM_ini, PM_deg)
% 1 seperator grid
PM=PM_ini;
PM.epsilonS(1)=PM.epsilonS_ini(1)+mean(PM_deg.delta_epsilonSN);
PM.epsilonE = PM.epsilonE_ini+PM_deg.delta_epsilonE;
PM.equiCapacity= PM.equiCapacity_ini*(1-PM_deg.AMloss); %*PM.factor_Q;
PM.RSEI = PM.RSEI_ini + max(PM_deg.delta_Rsei)*PM.sandwichArea;  % SEI resistance mOhm*cm^2 (0.378 mOhm * PM.sandwichArea) % EIS
PM.RCT = PM.RCT_ini + max(PM_deg.delta_RDL)*PM.sandwichArea; % Contact resistance % GA
% PM.maxCs=zeros(1,3);
% PM.maxCs(1)=PM.equiCapacity*3600/(PM.Stoi100_ini(1)-PM.Stoi0_ini(1))/PM.sandwichArea/PM.deltaN/PM.epsilonS_ini(1)/PM.F;
% PM.maxCs(3)=PM.equiCapacity*3600/(PM.Stoi0_ini(2)-PM.Stoi100_ini(2))/PM.sandwichArea/PM.deltaP/PM.epsilonS_ini(3)/PM.F;

%%
PM.stop=0;

PM.grids    = 25;
PM.gridsR   = 1;
PM.time=0;                  % Time, s

PM.L        = PM.deltaN+PM.deltaSep+PM.deltaP;
                            % Thickness of single layer, cm
PM.N        = [1,floor(PM.deltaN/PM.L*PM.grids+0.5),...
    floor((PM.deltaN+PM.deltaSep)/PM.L*PM.grids+0.5),PM.grids];


%%


PM.i02    = zeros(PM.grids,1);
PM.i02(PM.N(1):PM.N(2)) = PM.i0_k(1);
PM.i02(PM.N(3):PM.N(4)) = PM.i0_k(3);

PM.Rs2(PM.N(1):PM.N(2))=PM.Rs(1);
PM.Rs2(PM.N(3):PM.N(4))=PM.Rs(3);

                            % Grids Index
PM.deltaL=ones(PM.grids+1,1);
    PM.deltaL(PM.N(1))=0;
%     PM.deltaL(PM.N(1)+1:PM.N(2))=PM.deltaN/(PM.N(2)-PM.N(1));
%     PM.deltaL(PM.N(2)+1:PM.N(3))=PM.deltaSep/(PM.N(3)-PM.N(2));
%     PM.deltaL(PM.N(3)+1:PM.N(4))=PM.deltaP/(PM.N(4)-PM.N(3));
    PM.deltaL(PM.N(1):PM.N(2))=PM.deltaN/(PM.N(2)-PM.N(1)+1);
    PM.deltaL(PM.N(2)+1:PM.N(3)-1)=PM.deltaSep/(PM.N(3)-PM.N(2)-1);
    PM.deltaL(PM.N(3):PM.N(4))=PM.deltaP/(PM.N(4)-PM.N(3)+1);
    PM.deltaL(PM.N(4)+1)=0;
    
    PM.deltaLX = zeros(PM.grids,1);
for q = 2:PM.grids
    PM.deltaLX(q) = PM.deltaLX(q-1) + PM.deltaL(q);
end

                            % Thickness of N,Sep, P /grid, cm
PM.jLiArea=[0;ones(PM.N(2)-PM.N(1),1);zeros(PM.N(3)-PM.N(2)-1,1);ones(PM.N(4)-PM.N(3)+1,1);0];    

temp1       = PM.jLiArea.*PM.deltaL;
% PM.jLiDL    = (temp1(1:PM.grids)+temp1(2:PM.grids+1))/2;
PM.jLiDL    = temp1;
%% PROPERTY  
% Transfer coefficient
PM.alpha = [0.5 0 0.5];
% PM.alpha_a=0.5;             % Anodic charge transfer coefficient
% PM.alpha_c=0.5;             % Cathodic charge transfer coefficient

%% parameter f(temperature)
% PM.PM_Ds=PM.Ds; 
% PM.PM_RCT=PM.RCT; 
% PM.Ea_Ds=2.475e+004; %PM.Ds_0=5; %PM.Ds_0=1.5;
% PM.Ea_RCT=1.968e+004;
% PM.RCT=PM.PM_RCT*exp(PM.Ea_RCT/PM.R*(1/T_sim-1/PM.T0));
% PM.Ds(3)=PM.PM_Ds(3)*exp(-PM.Ea_Ds/PM.R*(1/T_sim-1/PM.T0));
% PM.Ds2=zeros(PM.grids,1); PM.Ds2(PM.N(1):PM.N(2))=PM.Ds(1); PM.Ds2(PM.N(3):PM.N(4))=PM.Ds(3);
% 
% % parameter
% PM.PM_Ds=PM.Ds; 
% PM.PM_De=PM.De;
% PM.PM_Rsei=PM.RSEI; 
% PM.Ea_Ds=4.475e+004; PM.Ds_0=1.5;
% PM.Ea_De=4.858e+004; PM.De_0=0.506;
% PM.Ea_R=3.796e+004;  PM.R_0=0.9474;
% 
% PM.Ds=PM.PM_Ds*PM.Ds_0*exp(-PM.Ea_Ds/PM.R*(1/T_sim-1/PM.T0));
% PM.Ds2=zeros(PM.grids,1); PM.Ds2(PM.N(1):PM.N(2))=PM.Ds(1); PM.Ds2(PM.N(3):PM.N(4))=PM.Ds(3);
% PM.De=PM.PM_De*PM.De_0*exp(-PM.Ea_De/PM.R*(1/T_sim-1/PM.T0));
% PM.RSEI=PM.PM_Rsei*PM.R_0*exp(PM.Ea_R/PM.R*(1/T_sim-1/PM.T0));


%% VARIABLE

PM.aS=3*PM.epsilonS./PM.Rs; % Specific surface area of electrode, m^-1
% PM.RSEI2=zeros(PM.N(2),1)
% PM.RSEI2(PM.N(1):PM.N(2),1) = PM.RSEI_ini + (PM_deg.delta_Rsei).*PM.factor_R;
% Array
    % Array epsilonE
PM.epsilonE2=ones(PM.grids+1,1);
    PM.epsilonE2(PM.N(1))=0;
    PM.epsilonE2(PM.N(1)+1:PM.N(2))=PM.epsilonE(1);
    PM.epsilonE2(PM.N(2)+1:PM.N(3)-1)=PM.epsilonE(2);
    PM.epsilonE2(PM.N(3):PM.N(4))=PM.epsilonE(3);
    PM.epsilonE2(PM.N(4)+1)=0; 
    % Array Ds
PM.Ds2=zeros(PM.grids,1);
    PM.Ds2(PM.N(1):PM.N(2))=PM.Ds(1).*(1+PM_deg.delta_epsilonSN/PM.epsilonS_ini(1));
    PM.Ds2(PM.N(3):PM.N(4))=PM.Ds(3);
    % Array effDe
PM.effDe=PM.De*PM.epsilonE2.^PM.porosityExponent;
    % Array effsigma
PM.effsigma2=ones(PM.grids+1,1);
    PM.effsigma2(PM.N(1))=0;
    PM.effsigma2(PM.N(1):PM.N(2))=PM.sigma(1)*PM.epsilonS_ini(1).*(1+PM_deg.delta_epsilonSN/PM.epsilonS_ini(1));
%     PM.effsigma2(PM.N(1):PM.N(2))=PM.sigma(1)*PM.epsilonS(1);
    PM.effsigma2(PM.N(2)+1:PM.N(3)-1)=PM.sigma(2)*PM.epsilonS(2);
    PM.effsigma2(PM.N(3):PM.N(4))=PM.sigma(3)*PM.epsilonS(3);
    PM.effsigma2(PM.N(4)+1)=0;
    % Array aS

PM.effsigma=PM.sigma.*PM.epsilonS;

PM.aS2=ones(PM.grids,1);
    PM.aS2(PM.N(1):PM.N(2))=PM.aS_ini(1).*(1+PM_deg.delta_epsilonSN(1:end)/PM.epsilonS_ini(1));
%         PM.aS2(PM.N(1):PM.N(2))=PM.aS(1);
    PM.aS2(PM.N(2)+1:PM.N(3)-1)=0;
    PM.aS2(PM.N(3):PM.N(4))=PM.aS(3);
PM.aS2(PM.N(1))=PM.aS2(PM.N(1)+1);

PM.Kappa = zeros(PM.grids,1);
tempCe = ones(PM.grids,1)*1.2e-3;
PM.Kappa(PM.N(1):PM.N(4)) = 15.8*tempCe.*exp(0.85*(1000*tempCe).^1.4);
% PM.Kappa(PM.N(1):PM.N(4))=15.8*tempCe.*exp(-13472*tempCe.^1.4); %MS
PM.effKappa = zeros(PM.grids,1);
PM.effKappa(PM.N(1):PM.N(2)) = PM.Kappa(PM.N(1):PM.N(2))*PM.epsilonE(1)^PM.porosityExponent;
PM.effKappa(PM.N(2)+1:PM.N(3)-1) = PM.Kappa(PM.N(2)+1:PM.N(3)-1)*PM.epsilonE(2)^PM.porosityExponent;
PM.effKappa(PM.N(3):PM.N(4)) = PM.Kappa(PM.N(3):PM.N(4))*PM.epsilonE(3)^PM.porosityExponent;
PM.effKappaD=2*PM.R*PM.T_sim*PM.effKappa/PM.F*(PM.tn-1);

%% Equilibrium potential
PM.Csrange_neg = PM.equiCapacity*3600/PM.sandwichArea/PM.deltaN/PM.epsilonS_ini(1)/PM.F; % Cs,ave(negative) range
PM.Csrange_pos = PM.equiCapacity*3600/PM.sandwichArea/PM.deltaP/PM.epsilonS_ini(3)/PM.F; % Cs,ave(positive) range

PM.startmidpoint = 0.5;
PM.startmidpoint_cat = 0.5;

% PM.startmidpoint = 0.60; % More aging
% PM.startmidpoint_cat = 0.65; % More aging

x_neg_0 = PM.startmidpoint - (0.5*PM.Csrange_neg)./PM.maxCs_ini(1);   % x_neg_SOC=0%
x_neg_100 = PM.startmidpoint + (0.5*PM.Csrange_neg)./PM.maxCs_ini(1); % x_neg_SOC=100%
x_pos_100 = PM.startmidpoint_cat - (0.5*PM.Csrange_pos)./PM.maxCs_ini(3); % x_pos_SOC=100%
x_pos_0 = PM.startmidpoint_cat + (0.5*PM.Csrange_pos)./PM.maxCs_ini(3);   % x_pos_SOC=0%

PM.Stoi100  = [x_neg_100, x_pos_100];
PM.Stoi0    = [x_neg_0, x_pos_0];

ocv_n=size(PM.ocv,1);       % Partition Index
x_neg_data=x_neg_0:(x_neg_100-x_neg_0)/(ocv_n-1):x_neg_100;
                            % Partition: (Cs_surf/Cs_max) of negative electrode
x_pos_data=x_pos_0:(x_pos_100-x_pos_0)/(ocv_n-1):x_pos_100;
                            % Partition: (Cs_surf/Cs_max) of positive electrode
    if isreal(x_neg_data)==0 || isreal(x_pos_data)==0
        PM.stop=1; disp('Imaginary stoi number')
    end
    
% f(x),df(x)/dx~x w/ extension
PM.LiMnO_data=zeros(3,ocv_n+2);
PM.LiMnO_data(1,2:ocv_n+1)=x_pos_data;
                        % (Cs_surf/Cs_max) of positive electrode 0~100%
                        % x
PM.LiMnO_data(1,1)=1;
PM.LiMnO_data(1,ocv_n+2)=0;
PM.LiMnO_data(2,2:ocv_n+1)=PM.ocv(ocv_n:-1:1)'+equi_LixC(x_neg_data);
                        % Equilibrium potiential of positive electrode=
                        % OCV+Equilibrium potiential of negative electrode 0~100%
                        % f(x)
PM.LiMnO_data(2,1)=PM.LiMnO_data(2,2)+((PM.LiMnO_data(2,3)-PM.LiMnO_data(2,2))/(PM.LiMnO_data(1,3)-PM.LiMnO_data(1,2)))*(PM.LiMnO_data(1,1)-PM.LiMnO_data(1,2));
PM.LiMnO_data(2,ocv_n+2)=PM.LiMnO_data(2,ocv_n)+((PM.LiMnO_data(2,ocv_n+1)-PM.LiMnO_data(2,ocv_n))/(PM.LiMnO_data(1,ocv_n+1)-PM.LiMnO_data(1,ocv_n)))*(PM.LiMnO_data(1,ocv_n+2)-PM.LiMnO_data(1,ocv_n));
% PM.LiMnO_data(2,1)=0;

PM.LiMnO_data(3,2:ocv_n+1)=(PM.LiMnO_data(2,3:ocv_n+2)-PM.LiMnO_data(2,2:ocv_n+1))./...
    (PM.LiMnO_data(1,3:ocv_n+2)-PM.LiMnO_data(1,2:ocv_n+1)); 
                        % df(x)/dx
PM.LiMnO_data(3,1)=PM.LiMnO_data(3,2);
PM.LiMnO_data(3,ocv_n+2)=PM.LiMnO_data(3,ocv_n+1);

%% effect of ion loss and active material loss XD
epsilonSN_temp = (PM.epsilonS(1)+PM.epsilonS(1)*(1-PM_deg.AMloss))/2; % indicates that ion loss and active material loss take place together
shiftX = ((PM_deg.ionloss)*3600)/(epsilonSN_temp*PM.deltaN*PM.sandwichArea*PM.maxCs(1)*PM.F);
% x=stoiYtoX(y,PM)-shiftX;y=y;
y_EOD= fzero(@(y) equi_LiMnO(y,PM)-equi_LixC(stoiYtoX(y,PM)-shiftX)-min(PM.ocv), PM.Stoi0(2) ); %+0.045 +0.2 50% +0.3 30%
x_EOD=stoiYtoX(y_EOD,PM)-shiftX;
y_EOC= fzero(@(y) equi_LiMnO(y,PM)-equi_LixC(stoiYtoX(y,PM)-shiftX)-max(PM.ocv),PM.Stoi100(2)); %-0.045
x_EOC=stoiYtoX(y_EOC,PM)-shiftX;
% PM.Stoi0    = [x_EOD y_EOD];
% PM.Stoi100  = [x_EOC y_EOC];
PM.Stoi0    = [x_EOD x_pos_0];
PM.Stoi100  = [x_EOC x_pos_100];

function [x]=stoiYtoX(y,PM)
x=-(y-PM.Stoi0(2))*(PM.epsilonS(3)*PM.deltaP*PM.maxCs(3))/(PM.epsilonS(1)*PM.deltaN*PM.maxCs(1))+PM.Stoi0(1);
end
%%
% aa=(ocv_n:-1:1)./100;
% figure; hold on;
% yyaxis right; ylabel('Ueq-/V'); 
% plot(equi_LixC(x_neg_data)); ylim([0 1])
% yyaxis left; ylabel('V'); 
% plot(PM.LiMnO_data(2,2:end-1))
% hold on; plot(PM.ocv(ocv_n:-1:1)','r-')
% plot(-equi_LixC(x_neg_data)+PM.LiMnO_data(2,2:end-1),'k--')
% legend('Ueq+','OCV','Ueq+-Ueq-'); 
% ylim([3 4.3]); xlim([0,100]); Plot_default; xlabel('SOC/%')
% 
% figure; hold on;
% nexttile; 
% plot(PM.LiMnO_data(1,2:end-1), PM.LiMnO_data(2,2:end-1));  xlim([0.1, 0.9]) 
% set(gca, 'XDir','reverse'); ylabel('Ueq+/V'); xlabel('Stoi_+');Plot_default
% 
% nexttile; 
% plot(x_neg_data, equi_LixC(x_neg_data)); xlim([0.1, 0.9]); 
% ylabel('Ueq-/V'); xlabel('Stoi_-');Plot_default
PM.PM_Ds=PM.Ds;
PM.PM_RCT=PM.RCT;
end