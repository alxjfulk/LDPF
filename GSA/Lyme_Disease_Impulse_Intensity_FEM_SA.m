function [SeF, SlF, IlF, SnF, InF, SaF, IaF, SmF, ImF] = Lyme_Disease_Impulse_Intensity_FEM_SA(LHSmatrix,j)

gif_yes = true;

addpath('../src_MMPDElab');

global piM miuM betaM piT sigmaT miuT betaT tauT gammaT miue K DM DA
global nuE nuL nuN nuA nuM

% dimensional parameters

    piM =  LHSmatrix(j,1);
    miuM = LHSmatrix(j,2);
    betaM = LHSmatrix(j,3);
    piT = LHSmatrix(j,4);
    sigmaT = LHSmatrix(j,5);
    miuT = LHSmatrix(j,6);
    betaT = LHSmatrix(j,7);
    tauT = LHSmatrix(j,8);
    gammaT = LHSmatrix(j,9);
    miue = LHSmatrix(j,10);
    K = LHSmatrix(j,11);
    DM = LHSmatrix(j,12);
    DA = LHSmatrix(j,13);
%Parameters for prescribed fire
    
    Pulses = LHSmatrix(j,14); %number of prescribed fire events
    tau = LHSmatrix(j,15); %time before next event
%Low intensity
%     nuE = 0;
%     nuA = 1 - 2/54;
%     nuN = 1 - 12/54;
%     nuL = 1 - 40/54;
%     nuM = 1 - 52/110;

% High intensity
  nuA = 1-(22- 12)/22;  
  nuN = 1-(3- 3)/3;  
  nuL = 1-(159- 118)/159;
  nuM = 1- 52/110; 
  nuE = 0;  

% compute initial meshes

   X = linspace(0,30,60)';
   tri = [(1:60-1)',(2:60)'];
   [Nv,d] = size(X);
   N = size(tri,1);

   Xi_ref = X;

% define nodes_fixed: corners are fixed

   tri_bf = [1;60];
   Nbf = length(tri_bf);
   nodes_fixed = [1;60];
   
% parameters for integration
   
   npde = 9;
   dt = 1e-5;
   mmpde_tau = 1e-3;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   M = repmat(reshape(eye(d,d),1,[]),Nv,1);
   t = 0.0;
   dtmax = 0.1;
      
% compute initial solutions
  
   U = u_initial(X);
   
   % for pde definition
   
   pdedef.bfMark = ones(Nbf,1);
   pdedef.bftype = zeros(Nbf,npde); % all Naumann BCs
   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;
   
% perform integration (MP)
   
   Se = zeros(5000,Nv);
   Sl = zeros(5000,Nv);
   Il = zeros(5000,Nv);
   Sn = zeros(5000,Nv);
   In = zeros(5000,Nv);
   Sa = zeros(5000,Nv);
   Ia = zeros(5000,Nv);
   Sm = zeros(5000,Nv);
   Im = zeros(5000,Nv);
   T = zeros(5000, 1);
  
   Se(1,:) = U(:,1)';
   Sl(1,:) = U(:,2)';
   Il(1,:) = U(:,3)';
   Sn(1,:) = U(:,4)';
   In(1,:) = U(:,5)';
   Sa(1,:) = U(:,6)';
   Ia(1,:) = U(:,7)';
   Sm(1,:) = U(:,8)';
   Im(1,:) = U(:,9)';

   DT = zeros(20000,4);
   
   tcpu = cputime;
   moving_mesh = false;
   n = 0;

%for loop here for impulse system
for i = 1:Pulses
    tf = t + tau;
    check = 1;
    dt = 1e-5;
    dt0 = dt;
while check == 1
   
      % move the mesh
      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([t,t+dt],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      else
         Xnew = X;
      end
      Xdot = (Xnew-X)/dt;
      
      % integrate physical PDEs
      [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef);
      
      % update
      n = n+1;
      X = X + dt0*Xdot;
      U = Unew;
      %if your solutions are going negative, uncommenting this can help avoid issues
      %related to estimation error:
      %U = max(U,0);
      t = t + dt0;
      dt = min(dtmax,dt1);
      NN = n+1;
      Se(n,:) = U(:,1)';
      Sl(n,:) = U(:,2)';
      Il(n,:) = U(:,3)';
      Sn(n,:) = U(:,4)';
      In(n,:) = U(:,5)';
      Sa(n,:) = U(:,6)';
      Ia(n,:) = U(:,7)';
      Sm(n,:) = U(:,8)';
      Im(n,:) = U(:,9)';
      T(n) = t;


      if (t+dt>tf), dt=tf-t; end
      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e',n,t,dt0,dt1);
      
      DT(n,:) = [t, dt0, norm(U(:,1),Inf), norm(U(:,2),Inf)];
      
      if (t>=tf-100*eps || dt < 100*eps), check = 0; end
      
end
 
%implement the pulse
   U(:,1) = Impulse1(X).*U(:,1);
   U(:,2) = Impulse3(X).*U(:,2);
   U(:,3) = Impulse3(X).*U(:,3);
   U(:,4) = Impulse4(X).*U(:,4);
   U(:,5) = Impulse4(X).*U(:,5);
   U(:,6) = Impulse5(X).*U(:,6);
   U(:,7) = Impulse5(X).*U(:,7);
   U(:,8) = Impulse6(X).*U(:,8);
   U(:,9) = Impulse6(X).*U(:,9);
end  
tcpu = cputime-tcpu;
fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
SeF = sum(Se(NN-1,:));
SlF = sum(Sl(NN-1,:));
IlF = sum(Il(NN-1,:));
SnF = sum(Sn(NN-1,:));
InF = sum(In(NN-1,:));
SaF = sum(Sa(NN-1,:));
IaF = sum(Ia(NN-1,:));
SmF = sum(Sm(NN-1,:));
ImF = sum(Im(NN-1,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PDE definition
function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    global piM miuM betaM piT sigmaT miuT betaT tauT gammaT miue K DM DA
    Nm = u(:,8) + u(:,9);
    lambdaT = betaT.*u(:,9)./(abs(Nm)+eps);
    lambdaM = betaM.*(u(:,3)+u(:,5)+u(:,7))./(abs(Nm)+eps);
    %pde model
    switch(ipde)
        case 1
        F = ut(:,1).*v(:) - piT.*(1-(u(:,1)/K)).*(u(:,6)+u(:,7)).*v(:) ...
            + (sigmaT + miue).*u(:,1).*v(:);
        case 2 %
        F = ut(:,2).*v(:) + DM.*du(:,2).*dv(:,1) - sigmaT.*u(:,1).*v(:) + lambdaT.*u(:,2).*v(:)...
            + (tauT+miuT).*u(:,2).*v(:);
        case 3 % 
        F = ut(:,3).*v(:) + DM.*du(:,3).*dv(:,1) - lambdaT.*u(:,2).*v(:) + (tauT+miuT).*u(:,3).*v(:);     
        case 4
        F = ut(:,4).*v(:) + DM.*du(:,4).*dv(:,1) - tauT.*u(:,2).*v(:) + lambdaT.*u(:,4).*v(:) ...
            + (gammaT+miuT).*u(:,4).*v(:);
        case 5 
        F = ut(:,5).*v(:) + DM.*du(:,5).*dv(:,1) - tauT.*u(:,3).*v(:) - lambdaT.*u(:,4).*v(:) ...
            + (gammaT+miuT).*u(:,5).*v(:);
        case 6 %
        F = ut(:,6).*v(:)+ DA.*du(:,6).*dv(:,1)  - gammaT.*u(:,4).*v(:) + lambdaT.*u(:,6).*v(:)...
            + miuT.*u(:,6).*v(:);
        case 7 %
        F = ut(:,7).*v(:)+ DA.*du(:,7).*dv(:,1)  - gammaT.*u(:,5).*v(:) - lambdaT.*u(:,6).*v(:)...
            + miuT.*u(:,7).*v(:);
        case 8 %
        F = ut(:,8).*v(:) + DM.*du(:,8).*dv(:,1) - piM.*v(:) + lambdaM.*u(:,8).*v(:) + miuM.*u(:,8).*v(:); 
        case 9 % 
        F = ut(:,9).*v(:) + DM.*du(:,9).*dv(:,1) - lambdaM.*u(:,8).*v(:) + miuM.*u(:,9).*v(:);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,ipde);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial conditions
function u = u_initial(X)

   u = zeros(size(X,1),9);
   %for uniform IC:
%    u(:,1) = 10000;
%    u(:,2) = 0;
%    u(:,3) = 100;
%    u(:,4) = 0;
%    u(:,5) = 5;
%    u(:,6) = 0;
%    u(:,7) = 100;
%    u(:,8) = 1000;
%    u(:,9) = 10;

%  for nonuniform IC:
   u(0 <= X(:,1) & X(:,1) <= 15,1) = 10000;
   u(:,2) = 0;
   u(0 <= X(:,1) & X(:,1) <= 15,3) = 100;
   u(:,4) = 0;
   u(0 <= X(:,1) & X(:,1) <= 15,5) = 500;
   u(:,6) = 0;
   u(0 <= X(:,1) & X(:,1) <= 15,7) = 100;
   u(:,8) = 1000;
   u(:,9) = 0;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fe = Impulse1(X)
global nuE
% X is in dimesnional form, Fl is a dimensionless proportion
   Fe = ones(size(X,1),1);
   Fe(:) = 1-nuE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Flm = Impulse3(X)
global nuL
% X is in dimesnional form, Flm is a dimensionless proportion
   Flm = ones(size(X,1),1);
   Flm(:) = 1-nuL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fn = Impulse4(X)
global nuN
% X is in dimesnional form, Ft is a dimensionless proportion
   Fn = ones(size(X,1),1);  
   Fn(:) = 1-nuN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fa = Impulse5(X)
global nuA

% X is in dimesnional form, Ft is a dimensionless proportion
   Fa = ones(size(X,1),1);  
   Fa(:) = 1-nuA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fm = Impulse6(X)
global nuM

% X is in dimesnional form, Ft is a dimensionless proportion
   Fm = ones(size(X,1),1);  
   Fm(:) = 1-nuM;