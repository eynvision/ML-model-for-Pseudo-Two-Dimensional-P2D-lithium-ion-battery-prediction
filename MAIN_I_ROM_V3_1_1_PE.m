function [Simdata,Simheat] =MAIN_I_ROM_V3_1_1_PE(Expdata,x,Crate)
load dOCVdT30KPA.mat
load DHDC.mat

    switch Crate
        case 1/3
            load GME101_p3C_EE_25oC_discharging.mat; 
        case 1/5
            load GME101_p2C_EE_25oC_discharging.mat;      
        case 1/1
            load GME101_1C_EE_25oC_discharging.mat;    
        case 1.5/1
            load GME101_1p5C_EE_25oC_discharging.mat;    
    end



format short g
global PM
%%  Input
%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%
Temperature = 25; % Ambient temperature C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn=length(Expdata.Time_S);

ambientT = Temperature+273.15; % Ambient temperature K
T_sim = ambientT; % Cell temperature K
DT_sim = Expdata.Time_S(10)-Expdata.Time_S(9); % Simulation Time Interval

PM=Parameters_sensitivity_pe(x);

Vt0 = Expdata.Voltage_V(1);
if Vt0 >= max(PM.ocv)
    SOC0 = 1;
elseif Vt0 <= min(PM.ocv)
    SOC0 = 0;
else
    SOC0 = interp1(PM.ocv, PM.soc,Vt0);
end

Cs_0 = [PM.startmidpoint.*PM.maxCs(1) - 0.5*PM.Csrange_neg,PM.startmidpoint.*PM.maxCs(3) + 0.5*PM.Csrange_pos]';
Cs0 = [SOC0.*PM.Csrange_neg+Cs_0(1),Cs_0(2)-SOC0.*PM.Csrange_pos]';




Cs_ave_n(PM.N(1):PM.N(2),1) = Cs0(1);
Cs_ave_p(1,1) = Cs0(2);

q_ave_n(PM.N(1):PM.N(2),1) = 0;
q_ave_p(1,1) = 0;

Cs_surf_n(PM.N(1):PM.N(2),1) = Cs0(1);
Cs_surf_p(1,1) = Cs0(2);

Ce = ones(PM.grids,1)*1.2e-3;
SOC = [SOC0,SOC0]; % Initial SOC
SOC_Ah = SOC0;

Phie = zeros(PM.grids,1);
Phis_n(PM.N(1):PM.N(2),1) = 0;
Phis_p = 0;

ref_stoi_n = Cs_ave_n(1,1)/PM.maxCs(1);
equiPotential(1) = equi_LixC(ref_stoi_n);
ref_stoi_p = Cs_ave_p(1,1)/PM.maxCs(3);
equiPotential(3) = equi_LiMnO(ref_stoi_p);
OCV = equiPotential(3) - equiPotential(1);

current = Expdata.Current_A(1);
jLi_n(PM.N(1):PM.N(2),1) = current/PM.sandwichArea/PM.deltaN;
jLi_p(1) =-current/PM.sandwichArea/PM.deltaP;

time=zeros(1,nn);
Vt=zeros(1,nn);
OCV=zeros(1,nn);
Heat_irr_prev=zeros(1,nn);
Heat_rev=zeros(1,nn);
Heat_tot=zeros(1,nn);
Heat_irr_anode_list=zeros(8,nn);
Heat_irr_sep_list=zeros(8,nn);
Heat_irr_cathode_list=zeros(8,nn);
Heat_irr_list=zeros(8,nn);

time(1) = 0;
Q(1) = 0;
Vt(1) = Vt0;
OCV(1) = Vt0;

%% Reduction of Ce
Ce0 = ones(PM.grids,1)*1.2e-3;
xCe = zeros(3,1);
[matrixA_Ce, matrixB_Ce, matrixC_Ce, matrixD_Ce] = i_matrix_Ce_RG;
ICe = eye(3);
matrixA_Ce_star = (ICe-matrixA_Ce/2*DT_sim)\(ICe+matrixA_Ce/2*DT_sim);
matrixB_Ce_star = (ICe-matrixA_Ce/2*DT_sim)\(matrixB_Ce*DT_sim);


%% Iteration
cnt = 1;
minVt=min(Expdata.Voltage_V);
% temp=find(Expdata.Voltage_V==minVt);
restflag=0;

while cnt <= length(Expdata.Time_S)-1
    cnt = cnt + 1;
    %     if cnt <= length(Expdata.Time_S)-1
    time(cnt) = Expdata.Time_S(cnt);
    DT_sim = time(cnt) - time(cnt-1);
    current = Expdata.Current_A(cnt);
    I(cnt) = current;
    jLi_n_ave(cnt) =  current/PM.sandwichArea/PM.deltaN;
    jLi_n(1,cnt) =  current/PM.sandwichArea/PM.deltaN;
    jLi_p(cnt) = -current/PM.sandwichArea/PM.deltaP;
    %% SOC Calculation
    Q(cnt)=Q(cnt-1)+current*DT_sim/3600;
    SOC_Ah(cnt)= SOC_Ah(cnt-1)-current*DT_sim/3600/PM.equiCapacity;
    %% Cs,ave
    k1_n = -3./PM.Rs_n./PM.F./PM.aS(1).*jLi_n(1,cnt);
    k1_p = -3./PM.Rs_p./PM.F./PM.aS(3).*jLi_p(1,cnt);
    Cs_ave_n(1,cnt) = Cs_ave_n(1,cnt-1)+k1_n.*DT_sim;
    Cs_ave_p(1,cnt) = Cs_ave_p(1,cnt-1)+k1_p.*DT_sim;
    %% q,ave
    q_ave_n(1,cnt) = ((1-15.*PM.Ds(1).*DT_sim./PM.Rs_n.^2).*q_ave_n(1,cnt-1)-45.*DT_sim./2./PM.Rs_n.^2./PM.F./PM.aS(1).*jLi_n(1,cnt))./(1+15.*PM.Ds(1).*DT_sim./PM.Rs_n.^2);
    q_ave_p(1,cnt) = ((1-15.*PM.Ds(3).*DT_sim./PM.Rs_p.^2).*q_ave_p(1,cnt-1)-45.*DT_sim./2./PM.Rs_p.^2./PM.F./PM.aS(3).*jLi_p(cnt))./(1+15.*PM.Ds(3).*DT_sim./PM.Rs_p.^2);
    %% Cs,surf
    Cs_surf_n(1,cnt) = Cs_ave_n(1,cnt) + 8.*PM.Rs_n./35.*q_ave_n(1,cnt) - PM.Rs_n./35./PM.Ds(1)./PM.F./PM.aS(1).*jLi_n(1,cnt);     Cs_surf_p(1,cnt) = Cs_ave_p(1,cnt) + 8.*PM.Rs_p./35.*q_ave_p(1,cnt) - PM.Rs_p./35./PM.Ds(3)./PM.F./PM.aS(3).*jLi_p(:,cnt);
    %% Calculation of electrolyte concentration
    xCe = matrixA_Ce_star*xCe+matrixB_Ce_star*current;
    Ce(:,cnt) = matrixC_Ce*xCe+matrixD_Ce*current+Ce0;
    %% ion conductivity
    PM.Kappa=zeros(PM.grids,1);
    tempCe=Ce(:,cnt);
    PM.Kappa(PM.N(1):PM.N(4)) = 15.8*tempCe.*exp(-13472*tempCe.^1.4);
    PM.effKappa = zeros(PM.grids,1);
    PM.effKappa((PM.N(1)):PM.N(2)) = PM.Kappa((PM.N(1)):PM.N(2))*PM.epsilonE(1)^PM.porosityExponent;
    PM.effKappa((PM.N(2)):PM.N(3)-1) = PM.Kappa((PM.N(2)):PM.N(3)-1)*PM.epsilonE(2)^PM.porosityExponent;
    PM.effKappa((PM.N(3)):PM.N(4)) = PM.Kappa((PM.N(3)):PM.N(4))*PM.epsilonE(3)^PM.porosityExponent;  
    PM.effKappaD=2*PM.R*T_sim*PM.effKappa/PM.F*(PM.tn-1);
    %% Phie
    Phie(1:PM.N(2),cnt) = -current.*PM.deltaLX(1:PM.N(2)).*PM.deltaLX(1:PM.N(2))./2./PM.effKappa(1:PM.N(2))./PM.deltaLX(PM.N(2))./PM.sandwichArea;
    Phie(PM.N(2):PM.N(3)-1,cnt) = -current.*PM.deltaLX(PM.N(2))./2./PM.effKappa(PM.N(2))./PM.sandwichArea -current.*(PM.deltaLX(PM.N(2):PM.N(3)-1)-PM.deltaLX(PM.N(2)))./PM.effKappa(PM.N(2):PM.N(3)-1)./PM.sandwichArea;
    Phie(PM.N(3):PM.N(4),cnt) = -current.*PM.deltaLX(PM.N(2))./2./PM.effKappa(PM.N(2))./PM.sandwichArea -current.*(PM.deltaLX(PM.N(3)-1)-PM.deltaLX(PM.N(2)))./PM.effKappa(PM.N(3)-1)./PM.sandwichArea -...
        current.*(PM.deltaLX(PM.N(3):PM.N(4))-PM.deltaLX(PM.N(3)-1)).*(2.*PM.deltaLX(PM.N(4))-PM.deltaLX(PM.N(3)-1)-PM.deltaLX(PM.N(3):PM.N(4)))./2./PM.effKappa(PM.N(3):PM.N(4))./(PM.deltaLX(PM.N(4))-PM.deltaLX(PM.N(3)-1))./PM.sandwichArea;
    %% Butler-Volmer Equation
    i0_n(1,1) = PM.k(1)*Ce(1,cnt).^PM.alpha(1)*(PM.maxCs(1)-Cs_surf_n(1,cnt)).^PM.alpha(1)*(Cs_surf_n(1,cnt)).^PM.alpha(3);
    i0_p = PM.k(3)*(mean(Ce(PM.N(3):PM.N(4),cnt))).^PM.alpha(1)*(PM.maxCs(3)-Cs_surf_p(cnt)).^PM.alpha(1)*(Cs_surf_p(cnt)).^PM.alpha(3);
    eta_n(1,cnt) = PM.R*T_sim./0.5./PM.F.*asinh(jLi_n(1,cnt)./2./PM.aS(1)./i0_n(1,1));
    eta_p(cnt) = PM.R*T_sim./0.5./PM.F.*asinh(jLi_p(cnt)./2./PM.aS(3)./i0_p);
    if isreal(eta_p(cnt))==0
        break;
    end
    
    %% OCV (Steady State)
    ref_stoi_n = Cs_ave_n(1,cnt)/PM.maxCs(1);
    equiPotential_n = equi_LixC(ref_stoi_n);
    ref_stoi_p = Cs_ave_p(cnt)/PM.maxCs(3);
    equiPotential_p = equi_LiMnO(ref_stoi_p);
    OCV(cnt) = equiPotential_p - equiPotential_n;
    %% Phis, anode (Transient State)
    stoi_n = Cs_surf_n(1,cnt)/PM.maxCs(1);
    stoi_p = Cs_surf_p(cnt)/PM.maxCs(3);
    Phis_n(1,cnt) = eta_n(1,cnt) + Phie(1,cnt) + equi_LixC(Cs_surf_n(1,cnt-1)/PM.maxCs(1)) + PM.RSEI*jLi_n(1,cnt)/PM.aS(1)*PM.sandwichArea;
    
    N = PM.N(2);
    dx = PM.deltaN./N;
    Phis_0 = Phis_n(1,cnt);
    
    Matrix_A = zeros(N+1,N);
    for i = 2:N
        Matrix_A(i+1,i) = -2;
    end
    for i = 2:N-1
        Matrix_A(i+1,i-1) = 1;
        Matrix_A(i+1,i+1) = 1;
    end
    Matrix_A(1,1) = 1;
    Matrix_A(2,2) = 2;
    Matrix_A(N+1,N-1) = 2;
    
    Matrix_B = zeros(N+1,1);
    Matrix_B(3:N,1) = dx*dx*jLi_n_ave(cnt)/PM.effsigma(1);
    Matrix_B(2,1) = 2*Phis_0+dx*dx*jLi_n_ave(cnt)/PM.effsigma(1)...
        - 2*dx*jLi_n_ave(cnt)/PM.effsigma(1)*PM.deltaN;
    Matrix_B(1,1) = Phis_0;
    
    Phis_n(:,cnt) = (Matrix_A'*Matrix_A)\Matrix_A'*Matrix_B;
    Phis_n(1,cnt) = Phis_0;
    
    %% Phis, cathode (Transient State)
    Phis_p(cnt) = eta_p(cnt) + mean(Phie(PM.N(3):PM.N(4),cnt)) + equi_LiMnO(stoi_p);
    %% Vt (Transit State)
    Vt(cnt) = Phis_p(cnt) -Phis_n(1,cnt) -PM.filmR*current;
    
    if isreal(Vt(cnt))==0
        break;
    end
    
%     if Vt(cnt)<2.4
%         break;
%     end
    
    %% State extention in Anode
    % Overpotential
    eta_n(1:PM.N(2),cnt) = Phis_n(1:PM.N(2),cnt) - Phie(1:PM.N(2),cnt) - equi_LixC(Cs_surf_n(1,cnt-1)./PM.maxCs(1)) - PM.RSEI*jLi_n(1,cnt)/PM.aS(1)*PM.sandwichArea;
    % Reaction current density
    i0_n(1:PM.N(2),1) = PM.k(1)*(Ce(1:PM.N(2),cnt)).^PM.alpha(1).*(PM.maxCs(1)-Cs_surf_n(1,cnt-1)).^PM.alpha(1).*(Cs_surf_n(1,cnt-1)).^PM.alpha(3);
    jLi_n(1:PM.N(2),cnt) = 2.*PM.aS(1).*i0_n(1:PM.N(2)).*sinh(0.5*PM.F./PM.R./T_sim.*eta_n(:,cnt));
    % Concentration
    % Cs,ave
    k1_n(1:PM.N(2)) = -3./PM.Rs_n./PM.F./PM.aS(1).*jLi_n(1:PM.N(2),cnt);
    Cs_ave_n(1:PM.N(2),cnt) = Cs_ave_n(1,cnt-1)+k1_n(1:PM.N(2))'.*DT_sim;
    % q,ave
    q_ave_n((1:PM.N(2)),cnt) = ((1-15.*PM.Ds(1).*DT_sim./PM.Rs_n.^2).*q_ave_n(1,cnt-1)-45.*DT_sim./2./PM.Rs_n.^2./PM.F./PM.aS(1).*jLi_n((1:PM.N(2)),cnt))./(1+15.*PM.Ds(1).*DT_sim./PM.Rs_n.^2);
    % Cs,surf
    Cs_surf_n(1:PM.N(2),cnt) = Cs_ave_n(1,cnt) + 8.*PM.Rs_n./35.*q_ave_n(1:PM.N(2),cnt) - PM.Rs_n./35./PM.Ds(1)./PM.F./PM.aS(1).*jLi_n(1:PM.N(2),cnt);
%     


    %% heat analysis
    jLi(PM.N(1):PM.N(2),1) = current/PM.sandwichArea/PM.deltaN;
    jLi(PM.N(3):PM.N(4),1) = -current/PM.sandwichArea/PM.deltaP;
    overPotential=[eta_n(:,cnt);0;eta_p(cnt)*ones(PM.N(4)-PM.N(3)+1,1)];
    Imicro1=PM.deltaL(PM.N(2))*sum(jLi(1:PM.N(2)));
    Imicro2=PM.deltaL(PM.N(3))*sum(jLi(PM.N(3):PM.N(4)));

    %% heat in electrode - electron transport
    phiS_out=[Phis_n(:,cnt);0;Phis_p(cnt)*ones(PM.N(4)-PM.N(3)+1,1)];
    Qirr_electrode=zeros(PM.grids,1);
    Qirr_electrode(PM.N(1))=Imicro1^2/PM.effsigma(1);

    Qirr_electrode(PM.N(1)+1:PM.N(2))=((phiS_out(PM.N(1)+1:PM.N(2))-phiS_out(PM.N(1):PM.N(2)-1))./PM.deltaL(PM.N(1)+1:PM.N(2))).^2*PM.effsigma(1);
    % only consider the anode because delta phi doesn't exist in cathode becauses of single particle

    lnCe=log(Ce(:,cnt));

    %% heat in electrolyte - ion transport (migration & diffusion)
    phiE_out= Phie(:,cnt);
    Qirr_electrolyte_cond=zeros(PM.grids,1);
    Qirr_electrolyte_diff=zeros(PM.grids,1);

    for i=2:PM.N(4)-1
        phiE1=(phiE_out(i)+phiE_out(i-1))/2;
        phiE2=(phiE_out(i)+phiE_out(i+1))/2;
        delL=(PM.deltaL(i)+PM.deltaL(i+1))/2;
        aveffK=(PM.effKappa(i)+PM.effKappa(i+1))/2;
        Qirr_electrolyte_cond(i)=((phiE2-phiE1)/delL).^2.*aveffK;%effKappa:ion  effective conductivity in electrolyte (migration)

        aveffKD=(PM.effKappaD(i)+PM.effKappaD(i+1))/2; %effKappaD: ion effectively diffusivity in electrolyte (diffusion)
        lnCe1=(lnCe(i)+lnCe(i-1))/2;
        lnCe2=(lnCe(i)+lnCe(i+1))/2;
        Qirr_electrolyte_diff(i)=aveffKD.*((lnCe2-lnCe1)./delL).*((phiE2-phiE1)./delL);
    end

    %% heat by activation potential
    Qirr_activation=zeros(PM.grids,1);%current overcome the overpotential
    Qirr_activation(PM.N(1):PM.N(4))=jLi(PM.N(1):PM.N(4)).*(overPotential(PM.N(1):PM.N(4)));%% heat by activation

    %% internal energy changed by the ion concentration causes heat
    Cs_surf=[Cs_surf_n(:,cnt);0;Cs_surf_p(cnt)*ones(PM.N(4)-PM.N(3)+1,1)];
    Cs_ave=[Cs_ave_n(:,cnt);0;Cs_ave_p(cnt)*ones(PM.N(4)-PM.N(3)+1,1)];
    Qirr_particle_anode2=(equi_LixC(Cs_surf(PM.N(1):PM.N(2))/PM.maxCs(1))-equi_LixC(Cs_ave(PM.N(1):PM.N(2))/PM.maxCs(1))).*jLi(PM.N(1):PM.N(2)); % enthalpy heating in anode %%My understanding: concetration diff within particle between the ave and surf → as current pass through the edge concentration is influenced → potential difference between the inter particle and surface particle → heat source
    Qirr_particle_cathode2=(equi_LiMnO(Cs_surf(PM.N(3):PM.N(4))/PM.maxCs(3))-equi_LiMnO(Cs_ave(PM.N(3):PM.N(4))/PM.maxCs(3))).*jLi(PM.N(3):PM.N(4));% enthalpy heating in cathod
    % https://www.sciencedirect.com/science/article/pii/S0306261921012320
    %annother title: Modeling and analysis of heat generation rate of a large format pouch-type lithium-ion battery considering degradation
    % As the electrons and ions meet at an active material site, chemical reactions take place at the surface of the active material particles, and the lithium ions are inserted into the particles or removed 
    %from the particles. As a result, a gradient in the concentration of the lithium ions at the surface of the particles is induced, which generates or consumes a large amount of power, dependent upon the 
    %direction of the reactions. Consequently, the chemical energy of the particles is changed, which is observed as changes in bond length between active material particles as lithium intercalates or deinteracalates 
    % and results in a change in the thermodynamic equilibrium of the active materials, resulting in a local change in the equilibrium potential of the active material. If the magnitude of the power consumed or generated 
    %during the change of the surface ion concentration and the magnitude of the change of the chemical energy is different, the remaining power is dissipated and produces additional heat

    %% heat in current collector
%     Qirr_CC=PM.filmR/PM.sandwichArea*current^2;
    Qirr_CC=PM.filmR*current^2;
%     Qirr_CC=19.4e-3*current^2;

    %% sum up
    PM.deltaL1=[PM.deltaL(PM.N(1)+1)/2;PM.deltaL(PM.N(1)+1:PM.N(4)-1);PM.deltaL(PM.N(4))/2];
    PM.deltaL1_electrode=[PM.deltaL(PM.N(1)+1)/2; PM.deltaL(PM.N(1)+1:PM.N(2)-1); PM.deltaL(PM.N(2))/2; 0;0;0; PM.deltaL(PM.N(3))/2 ; PM.deltaL(PM.N(3)+1:PM.N(4)-1); PM.deltaL(PM.N(4))/2];


    Qirr_electrode_anode=sum(Qirr_electrode(PM.N(1):PM.N(2)).*PM.deltaL1(PM.N(1):PM.N(2)));
    Qirr_electrode_cathode=sum(Qirr_electrode(PM.N(3)+1:PM.N(4)).*PM.deltaL1(PM.N(3)+1:PM.N(4)));

    Qirr_electrolyte_diff_anode=sum(Qirr_electrolyte_diff(PM.N(1):PM.N(2)).*PM.deltaL1(PM.N(1):PM.N(2)));
    Qirr_electrolyte_diff_sep=sum(Qirr_electrolyte_diff(PM.N(2)+1:PM.N(3)).*PM.deltaL1(PM.N(2)+1:PM.N(3)));
    Qirr_electrolyte_diff_cathode=sum(Qirr_electrolyte_diff(PM.N(3)+1:PM.N(4)).*PM.deltaL1(PM.N(3)+1:PM.N(4)));

    Qirr_electrolyte_cond_anode=sum(Qirr_electrolyte_cond(PM.N(1):PM.N(2)).*PM.deltaL1(PM.N(1):PM.N(2)));
    Qirr_electrolyte_cond_sep=sum(Qirr_electrolyte_cond(PM.N(2)+1:PM.N(3)).*PM.deltaL1(PM.N(2)+1:PM.N(3)));
    Qirr_electrolyte_cond_cathode=sum(Qirr_electrolyte_cond(PM.N(3)+1:PM.N(4)).*PM.deltaL1(PM.N(3)+1:PM.N(4)));

    Qirr_activation_anode=sum(Qirr_activation(PM.N(1):PM.N(2)).*PM.deltaL1(PM.N(1):PM.N(2)));
    Qirr_activation_cathode=sum(Qirr_activation(PM.N(3)+1:PM.N(4)).*PM.deltaL1(PM.N(3)+1:PM.N(4)));

    Qirr_particle_anode=sum(Qirr_particle_anode2.*PM.deltaL1(PM.N(1):PM.N(2))); % enthalpy heating in anode
    Qirr_particle_cathode=sum(Qirr_particle_cathode2.*PM.deltaL1(PM.N(3):PM.N(4)));
    % post calculation
    avgCs(1) = (sum(Cs_ave(PM.N(1):PM.N(2)))-(Cs_ave(PM.N(1))+Cs_ave(PM.N(2)))/2)/(PM.N(2)-PM.N(1));
    avgCs(2) = (sum(Cs_ave(PM.N(3):PM.N(4)))-(Cs_ave(PM.N(3))+Cs_ave(PM.N(4)))/2)/(PM.N(4)-PM.N(3));
    Stoi_neg = avgCs(1)/PM.maxCs(1);
    Stoi_pos = avgCs(2)/PM.maxCs(3);

    SOC_neg = (Stoi_neg-PM.Stoi_neg_0)./(PM.Stoi_neg_100-PM.Stoi_neg_0);
    %single particle hom
    Cs_surf_hom=[Cs_surf(1)*ones(PM.N(2),1);Cs_surf(PM.N(2)+1:PM.N(4))];
    Ctemp=(Cs_surf_hom-Cs_ave).^2;
    Cs_surf_prev=[Cs_surf_n(1,cnt-1)*ones(PM.N(2),1);0;Cs_surf_p(cnt-1)*ones(PM.N(4)-PM.N(3)+1,1)];
    Cs_ave_prev=[Cs_ave_n(:,cnt-1);0;Cs_ave_p(cnt-1)*ones(PM.N(4)-PM.N(3)+1,1)];
    Ctemp_prev=(Cs_surf_prev-Cs_ave_prev).^2;
        % from fitting
    % dUdT=ec_p(1).*SOC_neg .^7+ec_p(2).*SOC_neg .^6+ec_p(3).*SOC_neg .^5+ec_p(4).*SOC_neg .^4+ec_p(5).*SOC_neg .^3+ec_p(6).*SOC_neg .^2+ec_p(7).*SOC_neg+ec_p(8);
    % 
     dUdT= interp1(combined_data(:,2),combined_data(:,1),real(SOC_neg),'pchip');
    eccal2 = interp1(DHDC(:,2),DHDC(:,1),real(SOC_neg),'pchip');
    dHdC= eccal2;
%     dHdC=PM.F*(2267.*SOC_neg.^6-7555.*SOC_neg.^5+9873.*SOC_neg.^4-6393.*SOC_neg.^3+2124.*SOC_neg.^2-333.2.*SOC_neg.^1+18.93);
    HOM=abs(1/2*dHdC*(Ctemp-Ctemp_prev)/DT_sim);
    sumHOM=sum(HOM.*PM.deltaL1.*PM.sandwichArea);

    Qirr_mix_anode=sum(abs(1/2*dHdC*(Ctemp(1:PM.N(2))-Ctemp_prev(1:PM.N(2)))/DT_sim.*PM.deltaL1(1:PM.N(2)))); % heat of mixing in anode
    Qirr_mix_cathode=sum(abs(1/2*dHdC*(Ctemp(PM.N(3):PM.N(4))-Ctemp_prev(PM.N(3):PM.N(4)))/DT_sim.*PM.deltaL1(PM.N(3):PM.N(4))));% heat of mixing in cathod

    Heat_irr_anode=[Qirr_electrode_anode;Qirr_electrolyte_diff_anode;Qirr_electrolyte_cond_anode;Qirr_activation_anode;0;Qirr_particle_anode;Qirr_mix_anode;0]*PM.sandwichArea;
    Heat_irr_sep=[0;Qirr_electrolyte_diff_sep;Qirr_electrolyte_cond_sep;0;0;0;0;0]*PM.sandwichArea;
    Heat_irr_cathode=[Qirr_electrode_cathode;Qirr_electrolyte_diff_cathode;Qirr_electrolyte_cond_cathode;Qirr_activation_cathode;0;Qirr_particle_cathode;Qirr_mix_cathode;0]*PM.sandwichArea;

       



    Heat_SEI2=zeros(PM.N(2),1);
    Heat_SEI2(PM.N(1):PM.N(2))=PM.Rsei*jLi(PM.N(1):PM.N(2))./PM.aS2(PM.N(1):PM.N(2)).*jLi(PM.N(1):PM.N(2));

    Heat_SEI=sum(Heat_SEI2(PM.N(1):PM.N(2)).*PM.deltaL1(PM.N(1):PM.N(2))).*PM.sandwichArea ;


    Heat_irr=Heat_irr_anode+Heat_irr_sep+Heat_irr_cathode+[0;0;0;0;Qirr_CC;0;0;0]+[0;0;0;0;0;0;0;Heat_SEI];
    Heat_irr_anode_l=Heat_irr_anode+[0;0;0;0;Qirr_CC/2;0;0;0]+[0;0;0;0;0;0;0;Heat_SEI];
    Heat_irr_cathode_l=Heat_irr_cathode+[0;0;0;0;Qirr_CC/2;0;0;0];
%% lumped heat analysis
%     
    Heat_irr_prev(cnt) = current*(OCV(cnt)-Vt(cnt));
    
    Hear_rev_jli=zeros(PM.grids,1);
    Hear_rev_jli(PM.N(1):PM.N(2))=-jLi(PM.N(1):PM.N(2)).*T_sim*dUdT; 
    Hear_rev_jli(PM.N(3):PM.N(4))=jLi(PM.N(3):PM.N(4)).*T_sim*dUdT; 

    Heat_rev(cnt)=-current*T_sim*dUdT;
    Heat_tot(cnt)=Heat_rev(cnt)+sum(Heat_irr);
    
    Heat_irr_list(:,cnt)=Heat_irr;
    Heat_irr_anode_list(:,cnt)=Heat_irr_anode_l;
    Heat_irr_sep_list(:,cnt)=Heat_irr_sep;
    Heat_irr_cathode_list(:,cnt)=Heat_irr_cathode_l;


end
%% Save Data
Simdata.t = time';
Simdata.I = I';
Simdata.Vt = Vt';
Simdata.OCV = OCV';
Simdata.Heat_irr_anode=Heat_irr_anode_list;
Simdata.Heat_irr_sep=Heat_irr_sep_list;
Simdata.Heat_irr_cathode=Heat_irr_cathode_list;
Simdata.Heat_irr=Heat_irr_list;
Simdata.Heat_irr_tot=sum(Heat_irr_list);
Simdata.Heat_irr_lump=Heat_irr_prev;
Simdata.Heat_rev=Heat_rev;

Simdata.Heat_tot=Heat_tot;