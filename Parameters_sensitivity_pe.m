function [PM]= Parameters_sensitivity_pe(x,x_deg)
% x(1) = Dsn, x(2) = Dsp, x(3) = Rfilm
global PM


PM.deltaN   = 129.8e-4;  % negative thickness % provided by GM
PM.deltaSep   = 12e-4;   % Seperator thickness % provided by GM
PM.deltaP   = 85.9e-4;  % positive thickness % provided by GM
% 
% % 
% PM.epsilonS=[x(2), 0, x(1)];
% PM.epsilonE = [x(4),0.53,x(3)];
% PM.Rs= [x(6), 0, x(5)];
% PM.Ds=[x(8), 0, x(7)];
% PM.De=x(9);
% PM.filmR=x(10);
% x=[7.76638150982149e-11,1.16081618834729e-10,2.55646248309225e-07,5.08412129425597];
PM.epsilonS=[x(4), 0, x(8)];  % provided by GM
PM.epsilonE = [0.325,0.53,0.24]; % provided by GM
% PM.Rs= [9e-6, 0, 5e-6]; % provided by GM
PM.Rs= [x(5), 0, x(6)]; % provided by GM
PM.Ds=[x(1), 0, x(2)];
PM.De=x(3);
PM.filmR=x(9);

RSEI_EIS=x(10);
PM.RSEI = x(10);

PM.maxCs        = [0.033,0, 0.051];

PM.porosityExponent=1.5;
PM.tn=0.363;

PM.sigma=[1, 0, 0.1];

PM.i0_k  = [13.2,0,6.79];%M. Xiao, S.-Y. Choe / Journal of Power Sources 218 (2012) 357e367
PM.k= [13.2,0,6.79];


% %% Aging parameters
% PM.i0k_side=x_deg(1); % Exchange current density coefficient of SR
% PM.i0k_p=x_deg(2); %Exchange current density of LiPS
% 
% PM.Vsei = 100; % molar volume of SEI [cm3/mol]
% PM.Vp = 100; % molar volume of LiP [cm3/mol]
% PM.VDL = 7560; % molar volume of Deposit Layer [cm3/mol] by Literature [R.Fu]
% PM.Ve = 325; %%molar volume of electrolyte [cm3/mol] by Literature [R.Fu]
% 
% PM.alpha_c_side = .7; % Cathodic symmetric factor of SR
% PM.alpha_a_lp=0.33; % Anodic symmetric factor of LiP
% PM.alpha_c_lp=0.67; % Cathodic symmetric factor of LiP
% PM.equiPotential_side = 0.4; % Side reaction equilibrium potential [V]
% PM.equiPotential_lp = 0; % Lithium plating equilibrium potential [V]
% 
% PM.resistivitySEI = x_deg(4); %Ionic conductivity of SEI
% PM.resistivityDL = x_deg(5); %Ionic conductivity of DL
% PM.resistivityLP = 0;  %Ionic conductivity of LiP
% PM.lambda = 0.5; % Ratio for the formation of the secondary SEI from the plated Li
% 
% 




%% PARTITIPN
PM.grids=25;                % Grid number in cell
%% CONSTANT
PM.R=8.314;
PM.F=96485;
%% Geometry
PM.Length    = 54.5;                  % Length cm
PM.Width    = 11.1;                 % Width cm
PM.d    = 0.9;           % thickness [cm]
PM.V    = PM.Length*PM.Width*PM.d; % volume    [cm^3]

PM.layers_neg=40;
PM.layers_pos=38;
PM.layers   = PM.layers_neg+PM.layers_pos;

PM.sandwichArea = PM.Length*PM.Width*PM.layers;

% Mesh grids
PM.grids    = 25;
PM.gridsR   = 25;


PM.time=0;                  % Time, s
PM.deltat=0.5;              % Time interval (initial), s

PM.L        = PM.deltaN+PM.deltaSep+PM.deltaP;
                            % Thickness of single layer, cm
PM.N        = [1,floor(PM.deltaN/PM.L*PM.grids+0.5),...
    floor((PM.deltaN+PM.deltaSep)/PM.L*PM.grids+0.5),PM.grids];

PM.Rs2(PM.N(1):PM.N(2))=PM.Rs(1);
PM.Rs2(PM.N(3):PM.N(4))=PM.Rs(3);

PM.Rs_n=PM.Rs(1);
PM.Rs_p=PM.Rs(3);

                            % Grids Index
PM.deltaL=ones(PM.grids+1,1);
    PM.deltaL(PM.N(1))=0;
    PM.deltaL(PM.N(1)+1:PM.N(2))=PM.deltaN/(PM.N(2)-PM.N(1));
    PM.deltaL(PM.N(2)+1:PM.N(3))=PM.deltaSep/(PM.N(3)-PM.N(2));
    PM.deltaL(PM.N(3)+1:PM.N(4))=PM.deltaP/(PM.N(4)-PM.N(3));
    PM.deltaL(PM.N(4)+1)=0;
    
PM.deltaLX = zeros(PM.grids,1);
for q = 2:PM.grids+1
    PM.deltaLX(q) = PM.deltaLX(q-1) + PM.deltaL(q);
end
                            % Thickness of N,Sep, P /grid, cm
% PM.jLiArea=[0;ones(PM.N(2)-PM.N(1),1);zeros(PM.N(3)-PM.N(2),1);ones(PM.N(4)-PM.N(3),1);0];%%rom
                            % Calculation Weight
% Temperature
PM.T0=25+273.15;            % Temperature in ambient environment, K
%% PROPERTY  
% Transfer coefficient
PM.alpha_a=0.5;             % Anodic charge transfer coefficient
PM.alpha_c=0.5;             % Cathodic charge transfer coefficient
PM.alpha = [0.5 0 0.5];
%% VARIABLE

PM.aS=3*PM.epsilonS./PM.Rs; % Specific surface area of electrode, m^-1

% Array
    % Array epsilonE
PM.epsilonE2=ones(PM.grids+1,1);
    PM.epsilonE2(PM.N(1))=0;
    PM.epsilonE2(PM.N(1)+1:PM.N(2))=PM.epsilonE(1);
    PM.epsilonE2(PM.N(2)+1:PM.N(3))=PM.epsilonE(2);
    PM.epsilonE2(PM.N(3)+1:PM.N(4))=PM.epsilonE(3);
    PM.epsilonE2(PM.N(4)+1)=0; 
    % Array Ds
PM.Ds2=zeros(PM.grids,1);
    PM.Ds2(PM.N(1):PM.N(2))=PM.Ds(1);
    PM.Ds2(PM.N(3):PM.N(4))=PM.Ds(3);
    % Array effDe
PM.effDe=PM.De*PM.epsilonE2.^PM.porosityExponent;
    % Array effsigma
PM.effsigma=PM.sigma.*PM.epsilonS;

PM.aS2=ones(PM.grids,1);
    PM.aS2(PM.N(1):PM.N(2))=PM.aS(1);
    PM.aS2(PM.N(2)+1:PM.N(3)-1)=0;
    PM.aS2(PM.N(3):PM.N(4))=PM.aS(3);

%% Equilibrium potential
% if Crate==0.3
%     PM.equiCapacity=111.8317015731236; %1
%     PM.equiCapacity=109.9814843035368; %2
%     PM.equiCapacity=112.0408205815366; %3 using dch data
     PM.equiCapacity=x(7); % using ch data
% elseif Crate==0.25
%     PM.equiCapacity=111.7140450162599; %1
%     PM.equiCapacity= 109.9128558616463; %2
%     PM.equiCapacity=111.7371096364620; %3 using dch data
%     PM.equiCapacity=112.5787025838193; % using ch data
% end
    
load ocv_soc_GM_measured_ch.mat
PM.ocv=ocv_soc(:,1);               % OCV~SOC (100%~0%)
PM.soc=ocv_soc(:,2); 
                            
PM.Csrange_neg = PM.equiCapacity*3600/PM.sandwichArea/PM.deltaN/PM.epsilonS(1)/PM.F; % Cs,ave(negative) range
PM.Csrange_pos = PM.equiCapacity*3600/PM.sandwichArea/PM.deltaP/PM.epsilonS(3)/PM.F; % Cs,ave(positive) range

PM.startmidpoint = 0.50;
x_neg_0 = PM.startmidpoint - (0.5*PM.Csrange_neg)./PM.maxCs(1);   % x_neg_SOC=0%
x_neg_100 = PM.startmidpoint + (0.5*PM.Csrange_neg)./PM.maxCs(1); % x_neg_SOC=100%
x_pos_100 = PM.startmidpoint - (0.5*PM.Csrange_pos)./PM.maxCs(3); % x_pos_SOC=100%
x_pos_0 = PM.startmidpoint + (0.5*PM.Csrange_pos)./PM.maxCs(3);   % x_pos_SOC=0%


PM.x_neg_0 = x_neg_0;
PM.x_neg_100 = x_neg_100;
PM.x_pos_100 = x_pos_100;
PM.x_pos_0 = x_pos_0;
PM.Stoi_neg_100 = x_neg_100;
PM.Stoi_neg_0 = x_neg_0;
PM.Stoi_pos_100 = x_pos_100;
PM.Stoi_pos_0 = x_pos_0;

ocv_n=size(PM.ocv,1);       % Partition Index
x_neg_data=x_neg_0:(x_neg_100-x_neg_0)/(ocv_n-1):x_neg_100;
                            % Partition: (Cs_surf/Cs_max) of negative electrode
x_pos_data=x_pos_0:(x_pos_100-x_pos_0)/(ocv_n-1):x_pos_100;
                            % Partition: (Cs_surf/Cs_max) of positive electrode
                            
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

PM.LiMnO_data(3,2:ocv_n+1)=(PM.LiMnO_data(2,3:ocv_n+2)-PM.LiMnO_data(2,2:ocv_n+1))./...
    (PM.LiMnO_data(1,3:ocv_n+2)-PM.LiMnO_data(1,2:ocv_n+1)); 
                        % df(x)/dx
PM.LiMnO_data(3,2)=PM.LiMnO_data(3,2);                         
PM.LiMnO_data(3,1)=PM.LiMnO_data(3,2);
PM.LiMnO_data(3,ocv_n+2)=PM.LiMnO_data(3,ocv_n+1);

PM.Rsei = RSEI_EIS*PM.sandwichArea;                       
