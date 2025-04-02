function [out,Stop] = microMatrix_P2D_SPM(xCe,Cs_ave,q_ave,Cs_surf,current,jLi_tot,PM)
% global deltaTs
deltaT = PM.deltaTs;
T = PM.T_sim;  kk=PM.grids-PM.N(3);
if current>=0
    PM.i0k_side=PM.i0k_side*1e-3; 
%     PM.i0k_p=PM.i0k_p*1e-10;
%     PM.i0k_p=PM.i0k_p*1e-5;
end
%% [Ce Cs_ave q_ave phiSE Cs_surf]
    equiPotential=zeros(PM.N(3),1);
     if sum(Cs_surf(PM.N(1):PM.N(2))>PM.maxCs(1))>0
          disp('Cs-,surf over max');%Stop=1;
      elseif sum(Cs_surf(PM.N(1):PM.N(2))<0)>0
%           disp('Cs-,surf < 0');%Stop=1;
      elseif sum(Cs_surf(PM.N(3))>PM.maxCs(3))>0
          disp('Cs+,surf over max');%Stop=1;
     end
      
    Cs_surf0=Cs_surf+2e-6;
    num = 0; Stop=0;
%     current=PM.sandwichArea*CVinput;
%% first particle's reaction rate  
    jLi_n =  current/PM.sandwichArea/PM.deltaN;
    jLi_p = -current/PM.sandwichArea/PM.deltaP;
%% Reduction of Ce
    Ce0 = ones(PM.grids,1)*1.2e-3;
    [matrixA_Ce, matrixB_Ce, matrixC_Ce, matrixD_Ce] = i_matrix_Ce_RG(PM);
    ICe = eye(3);
    matrixA_Ce_star = (ICe-matrixA_Ce/2*deltaT)\(ICe+matrixA_Ce/2*deltaT);
    matrixB_Ce_star = (ICe-matrixA_Ce/2*deltaT)\(matrixB_Ce*deltaT);
    xCe = matrixA_Ce_star*xCe+matrixB_Ce_star*current;
    Ce = matrixC_Ce*xCe+matrixD_Ce*current+Ce0;
    Ce=max(Ce,0);
    %%
    PM.Kappa = zeros(PM.grids,1);
tempCe = ones(PM.grids,1)*1.2e-3;
% tempCe=abs((Ce(PM.N(1):PM.N(4)-1)+Ce(PM.N(1)+1:PM.N(4)))/2);
% tempCe=Ce; %easy to get error
PM.Kappa(PM.N(1):PM.N(4)) = 15.8*tempCe.*exp(0.85*(1000*tempCe).^1.4);
% PM.Kappa(PM.N(1):PM.N(4))=15.8*tempCe.*exp(-13472*tempCe.^1.4);


PM.effKappa = zeros(PM.grids,1);
PM.effKappa(PM.N(1):PM.N(2)) = PM.Kappa(PM.N(1):PM.N(2))*PM.epsilonE(1)^PM.porosityExponent;
PM.effKappa(PM.N(2)+1:PM.N(3)-1) = PM.Kappa(PM.N(2)+1:PM.N(3)-1)*PM.epsilonE(2)^PM.porosityExponent;
PM.effKappa(PM.N(3):PM.N(4)) = PM.Kappa(PM.N(3):PM.N(4))*PM.epsilonE(3)^PM.porosityExponent;
PM.effKappaD=2*PM.R*PM.T_sim*PM.effKappa/PM.F*(PM.tn-1);

%%
    k_sep = PM.N(3)-PM.N(2)-1;
    idxCs_ave = [0 (PM.N(1):PM.N(2)) zeros(1,k_sep) PM.N(3)-k_sep 0];
    idxq_ave = [0 PM.grids-k_sep-kk+(PM.N(1):PM.N(2)) zeros(1,k_sep) (PM.grids+PM.N(3)-kk-2*k_sep) 0];
    idxS=[0 PM.grids*2-k_sep*2-kk*2+(PM.N(1):PM.N(2)) zeros(1,k_sep) PM.grids*2-k_sep*3-2*kk+PM.N(3) 0];
    idxCs=[0 PM.grids*3-k_sep*3-kk*3+(PM.N(1):PM.N(2)) zeros(1,k_sep) PM.grids*3-k_sep*4-3*kk+PM.N(3) 0];
    idxE=[0 PM.grids*4-k_sep*4-4*kk+(1:PM.grids) 0];
    indexMax=idxE(PM.N(4)+1);
    temp_surf = 0;
%     while (max(abs(Cs_surf-Cs_surf0)) > 1e-6)
        matrixA = zeros(indexMax,indexMax);
        matrixB = zeros(indexMax,indexMax);
%         matrixC = zeros(indexMax,1);
        matrixD = zeros(indexMax,1);
        matrixE = zeros(1,indexMax);
        matrixE(1,idxS(end-1)) = 1; % calculate the terminal voltage
        matrixE(1,idxS(PM.N(1)+1)) = -1;
        matrixF = -PM.RCT;
        xx = zeros(PM.N(4),1);
        xxx = [Cs_ave(PM.N(1):PM.N(2));Cs_ave(PM.N(3));q_ave(PM.N(1):PM.N(2));q_ave(PM.N(3));...
            xx(PM.N(1):PM.N(2));xx(PM.N(3));Cs_surf(PM.N(1):PM.N(2));Cs_surf(PM.N(3));xx(PM.N(1):PM.N(4))];
        Cs_surf0=Cs_surf;
        i0=zeros(PM.grids,1);
%         i0(PM.N(1):PM.N(2))=PM.i0_k(1)*(abs(Ce(PM.N(1):PM.N(2)))'.^PM.alpha(1).*abs((PM.maxCs(1)-Cs_surf0(PM.N(1):PM.N(2),1)')).^PM.alpha(1).*Cs_surf0(PM.N(1):PM.N(2),1)').^PM.alpha(3);
        i0(PM.N(1):PM.N(2))=PM.i0_k(1)*(Ce(PM.N(1):PM.N(2)).*(PM.maxCs(1)-Cs_surf0(PM.N(1):PM.N(2),1)).*Cs_surf0(PM.N(1):PM.N(2),1)).^0.5;
        i0(PM.N(3):PM.N(4))=PM.i0_k(3)*(mean(Ce(PM.N(3):PM.N(4)))'.*(PM.maxCs(3)-Cs_surf0(PM.N(3),1)').*Cs_surf0(PM.N(3),1)').^0.5;
        i0(PM.N(1):PM.N(2)) = i0(PM.N(1):PM.N(2))./(1+PM.F/PM.R/T*PM.RSEI.*i0(PM.N(1):PM.N(2)));
        x=Cs_surf0(PM.N(1):PM.N(2),1)/PM.maxCs(1);
%         x=Cs_ave(PM.N(1):PM.N(2),1)/PM.maxCs(1);
        equiPotential(PM.N(1):PM.N(2))=equi_LixC(x);
        y=Cs_surf0(PM.N(3),1)/PM.maxCs(3);
        y=real(y);
        equiPotential(PM.N(3))=equi_LiMnO(y,PM);
        temp_equi(PM.N(1):PM.N(2),1) = d_equi_LixC(x)/PM.maxCs(1);
        temp_equi(PM.N(3),1) = d_equi_LiMnO(y,PM)/PM.maxCs(3);
        temp_equi(PM.N(2)+1:PM.N(3)-1)=0;
        temp_U = equiPotential-temp_equi.*Cs_surf0;
        %% Cs_ave
        for i=1:PM.N(3)
            if idxCs_ave(i+1)~=0
                temp = 3*deltaT/PM.Rs2(i)/PM.R/T*i0(i);
                matrixA(idxCs_ave(i+1),idxCs_ave(i+1))=1; % 1-12
                matrixA(idxCs_ave(i+1),idxE(i+1)) = -temp;
                matrixA(idxCs_ave(i+1),idxS(i+1)) = temp;% 25-36
                matrixA(idxCs_ave(i+1),idxCs(i+1)) = -temp*temp_equi(i); % 37-48
                matrixB(idxCs_ave(i+1),idxCs_ave(i+1))=1;
                matrixD(idxCs_ave(i+1),1)=temp*temp_U(i);
            end
            if i==PM.N(3)
%                 matrixD(idxCs_ave(i+1),1)=temp*temp_U(i);
                for kki=PM.N(3):PM.N(4)
                     matrixA(idxCs_ave(i+1),idxE(kki+1))=-temp/(PM.N(4)-PM.N(3)+1);
                end
            end
        end
%         temp_jli=temp;

        %% q_ave
        for i=1:PM.N(3)
            if idxq_ave(i+1)~=0
                temp = 45/2*deltaT/PM.Rs2(i)^2/PM.R/T*i0(i);
                matrixA(idxq_ave(i+1),idxq_ave(i+1))=1+deltaT*30*PM.Ds2(i)/PM.Rs2(i)^2/2; % 
                matrixA(idxq_ave(i+1),idxE(i+1)) = -temp;
                matrixA(idxq_ave(i+1),idxS(i+1)) = temp;
                matrixA(idxq_ave(i+1),idxCs(i+1)) = -temp*temp_equi(i);
                matrixB(idxq_ave(i+1),idxq_ave(i+1))=(1-deltaT*30*PM.Ds2(i)/PM.Rs2(i)^2/2);
                matrixD(idxq_ave(i+1),1)=temp*temp_U(i);
            end
            if i==PM.N(3)
%                 matrixD(idxq_ave(i+1),1)=temp*temp_U(i);
                for kki=PM.N(3):PM.N(4)
                     matrixA(idxq_ave(i+1),idxE(kki+1))=-temp/(PM.N(4)-PM.N(3)+1);
                end
            end
        end   
        %% phiS equation
        matrixA(idxS(2),idxS(2))=1;
        for i=PM.N(1)+1:PM.N(2)
            if idxS(i+1)~=0
                temp1=(PM.jLiArea(i)*PM.deltaL(i)+PM.jLiArea(i+1)*PM.deltaL(i+1))/2*PM.aS2(i)*i0(i)*PM.F/PM.R/T;
                matrixA(idxS(i+1),idxS(i+1))=-temp1;
                matrixA(idxS(i+1),idxE(i+1))=temp1;
                matrixA(idxS(i+1),idxCs(i+1))=temp1*temp_equi(i);
                matrixD(idxS(i+1),1)=-temp1*temp_U(i);
                if idxS(i)~=0
                    matrixA(idxS(i+1),idxS(i))=PM.effsigma2(i)/PM.deltaL(i);
                    matrixA(idxS(i+1),idxS(i+1))=matrixA(idxS(i+1),idxS(i+1))-PM.effsigma2(i)/PM.deltaL(i);
                end
                if idxS(i+2)~=0
                    matrixA(idxS(i+1),idxS(i+2))=PM.effsigma2(i+1)/PM.deltaL(i+1);
                    matrixA(idxS(i+1),idxS(i+1))=matrixA(idxS(i+1),idxS(i+1))-PM.effsigma2(i+1)/PM.deltaL(i+1);
                end
            end
        end
        %% phiS in cathod particle only 1 particle
         i=PM.N(3);
         matrixA(idxS(i+1),idxS(i+1))=1;
         matrixA(idxS(i+1),idxCs(i+1))=temp_equi(i,1)*-1;
%          matrixD(idxS(i+1),1)=equiPotential(PM.N(3))+PM.R*T/PM.aS2(i)/i0(i)/PM.F*jLi_p;
         matrixD(idxS(i+1),1)=temp_U(i)+PM.R*T/PM.aS2(i)/i0(i)/PM.F*jLi_p;
            for kki=PM.N(3):PM.N(4)
                 matrixA(idxS(i+1),idxE(kki+1))=-1/(PM.N(4)-PM.N(3)+1);
            end
        %% Cs_surf equation
        for i=1:PM.N(3)
            if idxCs(i+1)~=0
                temp1=i0(i)/PM.R/T;
                matrixA(idxCs(i+1),idxCs_ave(i+1))=-35*PM.Ds2(i)/PM.Rs2(i);
                matrixA(idxCs(i+1),idxq_ave(i+1))=-8*PM.Ds2(i);
                matrixA(idxCs(i+1),idxS(i+1))=temp1;
                matrixA(idxCs(i+1),idxE(i+1))=-temp1;
                matrixA(idxCs(i+1),idxCs(i+1))=-temp1*(temp_equi(i))+35*PM.Ds2(i)/PM.Rs2(i);
                matrixD(idxCs(i+1),1)=temp1*(temp_U(i));
            end
            if i==PM.N(3)
                for kki=PM.N(3):PM.N(4)
                     matrixA(idxCs(i+1),idxE(kki+1))=-temp1/(PM.N(4)-PM.N(3)+1);
                end
            end
        end
         %% phiE equation
         for i=1:PM.grids
            if idxE(i)~=0
%                 matrixA(idxE(i+1),idxCe(i))=PM.effKappaD(i)/PM.deltaL(i)/Ce(i);
%                 matrixA(idxE(i+1),idxCe(i+1))=matrixA(idxE(i+1),idxCe(i+1))-PM.effKappaD(i)/PM.deltaL(i)/Ce(i);
                matrixA(idxE(i+1),idxE(i))=PM.effKappa(i)/PM.deltaL(i);
                matrixA(idxE(i+1),idxE(i+1))=matrixA(idxE(i+1),idxE(i+1))-PM.effKappa(i)/PM.deltaL(i);
                matrixB(idxE(i+1),idxE(i+1))=matrixB(idxE(i+1),idxE(i+1));
            end
            if idxE(i+2)~=0
%                 matrixA(idxE(i+1),idxCe(i+2))=PM.effKappaD(i+1)/PM.deltaL(i+1)/Ce(i+1);
%                 matrixA(idxE(i+1),idxCe(i+1))=matrixA(idxE(i+1),idxCe(i+1))-PM.effKappaD(i+1)/PM.deltaL(i+1)/Ce(i+1);
                matrixA(idxE(i+1),idxE(i+2))=PM.effKappa(i+1)/PM.deltaL(i+1);
                matrixA(idxE(i+1),idxE(i+1))=matrixA(idxE(i+1),idxE(i+1))-PM.effKappa(i+1)/PM.deltaL(i+1);
                matrixB(idxE(i+1),idxE(i+1))=matrixB(idxE(i+1),idxE(i+1));
            end
            if i<=PM.N(2)
                if idxS(i+1)~=0
                    temp1=(PM.jLiArea(i)*PM.deltaL(i)+PM.jLiArea(i+1)*PM.deltaL(i+1))/2*PM.aS2(i)*i0(i)*PM.F/PM.R/T;
                    matrixA(idxE(i+1),idxS(i+1))=temp1;
                    matrixA(idxE(i+1),idxE(i+1))=matrixA(idxE(i+1),idxE(i+1))-temp1;
                    matrixA(idxE(i+1),idxCs(i+1))=-temp1*temp_equi(i);
                    matrixB(idxE(i+1),idxE(i+1))=matrixB(idxE(i+1),idxE(i+1));
                    matrixD(idxE(i+1),1) = temp1*temp_U(i);
                end
            end
            if  i>=PM.N(3)
                temp1=(PM.jLiArea(i)*PM.deltaL(i)+PM.jLiArea(i+1)*PM.deltaL(i+1))/2;
                matrixD(idxE(i+1),1) = -1*temp1*jLi_p;
            end
        end

%%
%         X = matrixA\(matrixB*xxx+matrixC*CVinput+matrixD);
          X = matrixA\(matrixB*xxx+matrixD);
          if sum(isnan(X))>0 || sum(isreal(X))==0;
              disp('X imag or Nan'); Stop=1; out=0; return
          end
         %%
        Cs_ave_out = [X((PM.N(1):PM.N(2))); zeros(k_sep,1); X(PM.N(3)-k_sep)];
        q_ave_out = [X(PM.grids-k_sep-kk+(PM.N(1):PM.N(2))); zeros(k_sep,1); X(PM.grids-k_sep*2-kk+PM.N(3))];
        phiE_out = X(PM.grids*4-k_sep*4-4*kk+(1:PM.grids));
        phiS_out = [X(PM.grids*2-k_sep*2-kk*2+(PM.N(1):PM.N(2))); zeros(k_sep,1); X(PM.grids*2-k_sep*3-kk*2+PM.N(3))];
%         phiS_out = [X(PM.grids*2-k_sep*2-kk*2+(PM.N(1):PM.N(2)))+PM.RSEI./PM.aS2(PM.N(1):PM.N(2)).*jLi_n; zeros(k_sep,1); X(PM.grids*2-k_sep*3-kk*2+PM.N(3))];

        Cs_surf_out = [X(PM.grids*3-k_sep*3-kk*3+(PM.N(1):PM.N(2))); zeros(k_sep,1); X(PM.grids*3-k_sep*4-kk*3+PM.N(3)) ];
        Cs_surf = Cs_surf_out;
      if sum(Cs_surf(PM.N(1):PM.N(2))>PM.maxCs(1))>0
          disp('Cs-,surf over max'); %Stop=1;  %Cs_surf(PM.N(1):PM.N(2))=min(Cs_surf(PM.N(1):PM.N(2)),PM.maxCs(1)); %Stop=1;
      elseif sum(Cs_surf(PM.N(1):PM.N(2))<0)>0
%           disp('Cs-,surf < 0'); %Stop=1; %Cs_surf=max(Cs_surf,0.0001); %
      elseif sum(Cs_surf(PM.N(3))>PM.maxCs(3))>0
          disp('Cs+,surf over max');%Stop=1;  %Cs_surf(PM.N(3))=min(Cs_surf(PM.N(3)),PM.maxCs(3)); % Stop=1;
      elseif sum(Ce<0)>0
          disp('Ce < 0');%Stop=1;
      end
        
%%
        phiE_out1=zeros(length(phiS_out),1);
        phiE_out1(PM.N(1):PM.N(2))=phiE_out(PM.N(1):PM.N(2));
        phiE_out1(PM.N(3))=mean(phiE_out(PM.N(3):PM.N(4)));
        x=Cs_surf(PM.N(1):PM.N(2),1)/PM.maxCs(1);
        equiPotential(PM.N(1):PM.N(2))=equi_LixC(x);
        y=Cs_surf(PM.N(3),1)/PM.maxCs(3);
        equiPotential(PM.N(3))=equi_LiMnO(real(y),PM);
        overPotential=phiS_out-phiE_out1-equiPotential;  
        overPotential(PM.N(1):PM.N(2))=overPotential(PM.N(1):PM.N(2))-PM.RSEI./PM.aS2(PM.N(1):PM.N(2)).*jLi_tot(PM.N(1):PM.N(2));    

        %%
        i01=zeros(length(phiS_out),1);
        i01(PM.N(1):PM.N(2))=i0(PM.N(1):PM.N(2));
        i01(PM.N(3))=mean(i0(PM.N(3):PM.N(4)));
        jLi_out=PM.F/PM.R/T*PM.aS2(PM.N(1):PM.N(3)).*i01.*(overPotential);
    %% Calculate SOC
    avgCs(1) = mean(Cs_ave_out(PM.N(1):PM.N(2)),1);
    avgCs(2) = Cs_ave_out(PM.N(3));
    Stoi=[avgCs(1)/PM.maxCs(1) avgCs(2)/PM.maxCs(3)];
    SOC=(Stoi-PM.Stoi0)./(PM.Stoi100-PM.Stoi0);
    OCP=[equi_LixC(Stoi(1)),equi_LiMnO(Stoi(2),PM)];
    OCV=equi_LiMnO(Stoi(2),PM)-equi_LixC(Stoi(1));
            %%
%         if CVmode==1
            Vt=matrixE*X+matrixF*current./PM.sandwichArea ; %- PM.RSEI*jLi_n/PM.aS(1);
%         elseif CVmode==2
%             CVoutput=sum((PM.deltaL(PM.N(1):PM.N(2)).*PM.jLiArea(PM.N(1):PM.N(2))+PM.deltaL(PM.N(1)+1:PM.N(2)+1).*PM.jLiArea(PM.N(1)+1:PM.N(2)+1))/2.*jLi_out(PM.N(1):PM.N(2)));
%             CVoutput=CVoutput*PM.sandwichArea;
%         end
        
    %% heat analysis
    Heat_irr = current*(OCV-Vt);
    dUdT=interp1(PM.dU_dT(:,2),PM.dU_dT(:,1),SOC(1),'linear','extrap');
    Heat_rev=-current*PM.T_sim*dUdT; % W/battery
    Heat_tot=Heat_rev+Heat_irr;
    
    out.Heat_irr=Heat_irr;
    out.Heat_rev=Heat_rev;
    out.Heat_tot=Heat_tot;

          %% side reaction 
        overPotential_side(PM.N(1):PM.N(2),1)=phiS_out(PM.N(1):PM.N(2))-phiE_out(PM.N(1):PM.N(2))-PM.equiPotential_side; % ...
                        -PM.RSEI./PM.aS2(PM.N(1):PM.N(2)).*jLi_tot(PM.N(1):PM.N(2));
        overPotential_LP(PM.N(1):PM.N(2),1)=phiS_out(PM.N(1):PM.N(2))-phiE_out(PM.N(1):PM.N(2)); %...
                       -PM.RSEI./PM.aS2(PM.N(1):PM.N(2)).*jLi_tot(PM.N(1):PM.N(2));

        CE=1e-5*ones((PM.N(2)-PM.N(2)+1),1);
        i0_side=PM.i0k_side*sqrt(Cs_surf0(PM.N(1):PM.N(2)).*CE);
        jLi_side_out(PM.N(1):PM.N(2),1)=-PM.aS2(PM.N(1):PM.N(2)).*i0_side.*exp(-PM.alpha_c_side*2*PM.F/PM.R/T*overPotential_side(PM.N(1):PM.N(2)));
        j_lp_net(PM.N(1):PM.N(2),1)=PM.aS2(PM.N(1):PM.N(2)).*PM.i0k_p.*(exp(PM.alpha_a_lp*PM.F/PM.R/T*overPotential_LP(PM.N(1):PM.N(2)))-exp(-PM.alpha_c_lp*PM.F/PM.R/T*overPotential_LP(PM.N(1):PM.N(2))));

        %% LiS
if current<=0
    jLi_ps_out=min(j_lp_net,0);
    jLi_p_out=jLi_ps_out;
    jLi_s_out=zeros(PM.N(2),1);
elseif current>0
    jLi_ps_out=max(j_lp_net,0);
    jLi_s_out=jLi_ps_out;
    jLi_p_out=zeros(PM.N(2),1);
end

jLi_tot=jLi_out(PM.N(1):PM.N(2))+jLi_side_out+jLi_p_out;
        %%
out.xCe=xCe;
out.Cs_ave = Cs_ave_out;
out.q_ave= q_ave_out;
out.Cs_surf= Cs_surf; %at CC    
out.Vt=Vt;
out.Ce=Ce(:,1);
out.Phis=phiS_out(:,1);
out.Phie=phiE_out(:,1);
out.jLi=jLi_out(:,1);
out.jLi_tot=jLi_tot(:,1);
out.eta=overPotential(:,1);
out.Ueq_ave=OCP;
out.Ueq_surf=equiPotential(:,1);
out.stoi_ave=Stoi;
% out.stoi_surf=[stoi_n(1); stoi_p]; %at CC
out.SOC=SOC;
out.OCV=OCV;
out.i0=i0;
out.overPotential_side=overPotential_side;
out.j_side=jLi_side_out;
out.overPotential_lp=overPotential_LP;
out.jLi_p=jLi_p_out;
out.jLi_s=jLi_s_out;
out.jLi_ps=jLi_ps_out;

      
        %% iteration
%         if Stop==0
%             num = num+1;
%             if abs(max(abs(Cs_surf-Cs_surf0))-temp_surf)<2e-7
%                 Cs_surf = Cs_surf0;
%             end
%             temp_surf = max(abs(Cs_surf-Cs_surf0));
%         else
%             break;
%         end
%      end
end
