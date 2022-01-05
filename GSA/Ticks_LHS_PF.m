clc

clf
clear
%% Sample size N
runs=1000;

% %% LHS MATRIX  %

Ticks_Parameter_Settings_LHS;  

piM = LHS_Call((1-0.2)*piM_base, piM_base,(1+0.2)*piM_base, 0 ,runs,'unif');  %  
miuM  = LHS_Call((1-0.2)*miuM_base, miuM_base,(1+0.2)*miuM_base, 0 ,runs,'unif');  %  
betaM  = LHS_Call((1-0.2)*betaM_base, betaM_base,(1+0.2)*betaM_base, 0 ,runs,'unif');  %  
piT  = LHS_Call((1-0.2)*piT_base, piT_base,(1+0.2)*piT_base, 0 ,runs,'unif');  %  
sigmaT  = LHS_Call((1-0.2)*sigmaT_base, sigmaT_base,(1+0.2)*sigmaT_base, 0 ,runs,'unif');  %  
miuT  = LHS_Call((1-0.2)*miuT_base, miuT_base,(1+0.2)*miuT_base, 0 ,runs,'unif');  %  
betaT  = LHS_Call((1-0.2)*betaT_base, betaT_base,(1+0.2)*betaT_base, 0 ,runs,'unif');  %  
tauT  = LHS_Call((1-0.2)*tauT_base, tauT_base,(1+0.2)*tauT_base, 0 ,runs,'unif');  %  
gammaT  = LHS_Call((1-0.2)*gammaT_base, gammaT_base,(1+0.2)*gammaT_base, 0 ,runs,'unif');  %  
miue  = LHS_Call((1-0.2)*miue_base, miue_base,(1+0.2)*miue_base, 0 ,runs,'unif');  %  
K  = LHS_Call((1-0.2)*K_base, K_base,(1+0.2)*K_base, 0 ,runs,'unif');
DM  = LHS_Call((1-0.2)*DM_base, DM_base,(1+0.2)*DM_base, 0 ,runs,'unif');
DA  = LHS_Call((1-0.2)*DA_base, DA_base,(1+0.2)*DA_base, 0 ,runs,'unif');% 
Pulses = LHS_Call(0.1, Pulses_base, 0.5, 0 ,runs,'po10');
tau = LHS_Call(0.1, tau_base, 0.5, 0 ,runs,'poi1');     


%% LHS MATRIX and PARAMETER LABELS

LHSmatrix=[piM miuM betaM piT sigmaT miuT betaT tauT gammaT miue K DM DA Pulses tau];    

for j=1:runs %Run solution x times choosing different values    
      % here we call a modified version of the code used in other
      % scenarios. You must open Lyme_Disease_Impulse_Intensity_FEM_SA in
      % order to edit the final time, initial conditions, diffusion
      % coefficients, etc. Currently, it is set to a heterogeneous domain
      % with homogeneous initial conditions and high intenisty burns.
      [SeF, SlF, IlF, SnF, InF, SaF, IaF, SmF, ImF] = Lyme_Disease_Impulse_Intensity_FEM_SA(LHSmatrix,j); 
    
     Se_lhs(:,j)=SeF;
     Sl_lhs(:,j)=SlF;
     Il_lhs(:,j)=IlF;
     Sn_lhs(:,j)=SnF;  
     In_lhs(:,j)=InF; 
     Sa_lhs(:,j)=SaF;
     Ia_lhs(:,j)=IaF; 
     Sm_lhs(:,j)=SmF; 
     Im_lhs(:,j)=ImF; 
     
     sumIT(:,j)=IlF+InF+IaF;
     
end

% CALCULATE PRCC
alpha = 0.05;


%%%% Expression for SA, SL, SN, SM
Se = K.*(gammaT.*sigmaT.*tauT.*(piT - miuT)-miuT.*(miue.*(gammaT + miuT).*(tauT + miuT)  + miuT.*sigmaT.*(gammaT + miuT + tauT)))./(gammaT.*piT.*sigmaT.*tauT); 
Sl = K.*(gammaT.*sigmaT.*tauT.*(piT - miuT)-miuT.*(miuT.*sigmaT.*(gammaT + miuT + tauT) + miue.*(gammaT + miuT).*(tauT + miuT)))./(tauT.*(tauT + miuT).*gammaT.*piT); 
Sn = K.*(gammaT.*sigmaT.*tauT.*(piT - miuT)-miuT.*(miue.*(gammaT + miuT).*(tauT + miuT) + miuT.*sigmaT.*(gammaT + miuT + tauT)))./(gammaT.*piT.*(gammaT.*miuT + gammaT.*tauT + miuT.^2 + miuT.*tauT)); 
Sa = K.*(gammaT.*sigmaT.*tauT.*(piT - miuT)-miuT.*(miue.*(gammaT + miuT).*(tauT + miuT) + miuT.*sigmaT.*(gammaT + miuT + tauT)))./(piT.*(gammaT.*miuT + gammaT.*tauT + miuT.^2 + miuT.*tauT).*miuT); 
Sm = piM./miuM;

[prcc_In sign_In sign_In_label]=PRCC(LHSmatrix,In_lhs,1:length(time_points),PRCC_var,alpha);
%[prcc_IT sign_IT sign_IT_label]=PRCC(LHSmatrix,sumIT,1:length(time_points),PRCC_var,alpha);

%For saving the vector containing the total number of either infectious
%nymphs or infectious ticks at low or high intensity
% save('In_LI_heter.mat', 'In_lhs')
save('IT_HI_heter.mat', 'sumIT')
disp("PRCC Values_heter")
disp(prcc_In)
disp("sign In_heter")
disp(sign_In)
disp("sign_In_label_heter")
disp(sign_In_label)
%disp(prcc_IT)
