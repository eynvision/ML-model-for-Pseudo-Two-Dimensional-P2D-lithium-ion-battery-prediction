function [matrixA_star, matrixB_star, matrixC_star, matrixD_star] = i_matrix_Ce_RG
global PM
%% Define matrix A, B, C, D
matrixA = zeros(PM.grids,PM.grids);
matrixB = zeros(PM.grids,1);
for i = 1:PM.grids
    tempL = PM.effDe(i)/PM.deltaL(i)/((PM.deltaL(i)*PM.epsilonE2(i)+PM.deltaL(i+1)*PM.epsilonE2(i+1))/2);
    tempR = PM.effDe(i+1)/PM.deltaL(i+1)/((PM.deltaL(i)*PM.epsilonE2(i)+PM.deltaL(i+1)*PM.epsilonE2(i+1))/2);
    % uniform current density
    if i <= PM.N(2)
        tempC = (1-0.363)/PM.F/PM.sandwichArea/PM.deltaN/((PM.epsilonE2(i)+PM.epsilonE2(i+1))/2);
    elseif i > PM.N(3)
        tempC = -(1-0.363)/PM.F/PM.sandwichArea/PM.deltaP/((PM.epsilonE2(i)+PM.epsilonE2(i+1))/2);
    else
        tempC = 0;
    end
    matrixA(i,i) = -tempL-tempR;
    matrixB(i) = tempC;
    if i ~= 1
        matrixA(i,i-1) = tempL;
    end
    if i ~= PM.grids
        matrixA(i,i+1) = tempR;
    end
    if i == 1
        matrixA(i,i+1) = tempR;
        matrixA(i,i) = -tempR;
    end
    if i == PM.grids
        matrixA(i,i-1) = tempL;
        matrixA(i,i) = -tempL;
    end
end
matrixC = eye(PM.grids);
matrixD = zeros(PM.grids,1);
mid = round(PM.grids/2);  % middle of separator
matrixA_mid = matrixA(mid,:);
matrixB_mid = matrixB(mid,1);
for i=1:PM.grids
    matrixA(i,:) = matrixA(i,:)-matrixA_mid;
    matrixB(i,1) = matrixB(i,1)-matrixB_mid;
end
i = mid;
matrixA(i,i) = -1;
Z_ss = -matrixC*(matrixA\matrixB)+matrixD;

%% Define matrix A_hat, B_hat, C_hat, D_hat
[eig_Vec_R,eig_Val] = eig(matrixA);
Eig_values = diag(eig_Val);
[Eig_values,Index] = sort(Eig_values,'descend');
eig_Vec_R = eig_Vec_R(:,Index);  % right EigenVector
eig_Vec_L = inv(eig_Vec_R);  % left EigenVector
Residue = zeros(PM.grids,PM.grids);
for i = 1:PM.grids
    Residue(:,i) = matrixC*eig_Vec_R(:,i)*eig_Vec_L(i,:)*matrixB/Eig_values(i);
end
% matrixA_hat = diag(Eig_values);
% matrixB_hat = ones(PM.grids,1);
matrixD_hat = Z_ss;
for i = 1:PM.grids
    %     matrixC_hat(:,i) = Residue(:,i)*Eig_values(i);
    matrixD_hat = matrixD_hat+Residue(:,i);
end

%% Define matrix A_star, B_star, C_star, D_star
R_grouped = zeros(PM.grids,3);
Eig_grouped = zeros(PM.grids,3);
Eig_grouped_avr = zeros(3,1);
for i=1:PM.grids
    if Eig_values(i) > -1e+0
        R_grouped(:,1) = R_grouped(:,1)+Residue(:,i);
        Eig_grouped(:,1) = Eig_grouped(:,1)+Eig_values(i)*Residue(:,i);
    elseif Eig_values(i) > -1e+1
        R_grouped(:,2) = R_grouped(:,2)+Residue(:,i);
        Eig_grouped(:,2) = Eig_grouped(:,2)+ Eig_values(i)*Residue(:,i);
    else
        R_grouped(:,3)=R_grouped(:,3)+Residue(:,i);
        Eig_grouped(:,3)=Eig_grouped(:,3)+ Eig_values(i)*Residue(:,i);
    end
end
for i=1:3
    if rank(R_grouped(:,i))==1
        Eig_grouped_avr(i,1) = Eig_grouped(:,i)'/R_grouped(:,i)';
    else
        Eig_grouped_avr(i,1) = 0;
    end
end
matrixA_star = diag(Eig_grouped_avr);
matrixB_star = ones(3,1);
matrixC_star = zeros(PM.grids,3);
for i=1:3
    matrixC_star(:,i) = R_grouped(:,i)*Eig_grouped_avr(i);
end
matrixD_star = matrixD_hat;
end