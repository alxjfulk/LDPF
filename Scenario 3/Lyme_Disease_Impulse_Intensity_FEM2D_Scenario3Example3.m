function Lyme_Disease_Impulse_Intensity_FEM2D_Scenario3Example3(jmax, kmax)

gif_yes = true;

addpath('../src_MMPDElab');

global piM miuM betaM piT sigmaT miuT betaT tauT gammaT miue K
global nuE nuL nuN nuA nuM x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2

% dimensional parameters

    piM = 0.02;
    miuM = 0.01;
    betaM = 0.9;
    piT = 456.36;
    sigmaT = 0.06670;
    miuT = 0.015;
    betaT = 0.9;
    tauT = 0.06180;
    gammaT = 0.04910;
    miue = 0.0025;
    K = 5000;

%Parameters for prescribed fire
    
    Pulses = 10; %number of prescribed fire events
    tau = 1; %time before next event
%Low intensity
    nuE = 0;
%     nuA = 1 - 2/54;
%     nuN = 1 - 12/54;
%     nuL = 1 - 40/54;
%     nuM = 1 - 52/110;

% % High intensity
  nuA = 1-(22- 12)/22;  
  nuN = 1-(3- 3)/3;  
  nuL = 1-(159- 118)/159;
  nuM = 1- 52/110; 


% compute initial meshes

   [X,tri] = MovMesh_rect2tri(linspace(0,30,jmax),linspace(0,30,kmax),1);
   TR = triangulation(tri,X);
   [Nv,d] = size(X);
   N = size(tri,1);
   Xi_ref = X;
   % define nodes_fixed: corners are fixed
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);
   nodes_fixed = unique(tri_bf);

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
   Se = zeros(10000,Nv);
   Sl = zeros(10000,Nv);
   Il = zeros(10000,Nv);
   Sn = zeros(10000,Nv);
   In = zeros(10000,Nv);
   Sa = zeros(10000,Nv);
   Ia = zeros(10000,Nv);
   Sm = zeros(10000,Nv);
   Im = zeros(10000,Nv);
   T = zeros(10000, 1);
  
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
   h = 1;
   
   tcpu = cputime;
   moving_mesh = false;
n = 0;

%implement a for loop here for impulse system
for i = 1:Pulses
    if i >= 10
        tau = 11;
    end
    tf = t + tau;
    check = 1;
    dt = 1e-5;
    dt0 = dt;
    x1loc1 = 0;
    x1loc2 = 10;
    y1loc1 = 0;
    y1loc2 = 10;
    x2loc1 = 0;
    x2loc2 = 10;
    y2loc1 = 0;
    y2loc2 = 10;
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
      U = max(U,0);
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

      figure(1)
      clf
      trisurf(tri,X(:,1),X(:,2),U(:,5))
      shading interp
      title(['Infectious Nymphs, t =' num2str(t)],'FontSize', 24)
      set(gca, 'FontSize', 20)
      xlabel('x')
      ylabel('y')
      view(0,90)
      %colormap(cm)
      zzaxis = colorbar;
      zzaxis.Label.String = 'Number of Ticks';
      zzaxis.Label.Rotation = -90;
      zzaxis.Label.Position(1) = 4.7;
      zzaxis.Label.Position(2) = 500;
      caxis([0,1000])
      zzaxis.LineWidth = 0.5;
      drawnow
      if gif_yes
      basename = 'INanimated11';
      filename = [num2str(h),basename,'.gif'];
      frame = getframe(gcf);
      im = frame2im(frame);
      [INgif, cm] = rgb2ind(im,256);
      if n == 1
          imwrite(INgif, cm , filename, 'gif', 'Loopcount', inf, 'DelayTime', .1);
      else
          imwrite(INgif, cm , filename, 'gif', 'WriteMode', 'append', 'DelayTime', .1);
      end
      end

      if (t+dt>tf), dt=tf-t; end
      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e',n,t,dt0,dt1);
      
      DT(n,:) = [t, dt0, norm(U(:,1),Inf), norm(U(:,2),Inf)];
      
      if (t>=tf-100*eps || dt < 100*eps), check = 0; end
      
end
 
%implement the pulse
   U(:,1) = Impulse1(X,t).*U(:,1);
   U(:,2) = Impulse3(X,t).*U(:,2);
   U(:,3) = Impulse3(X,t).*U(:,3);
   U(:,4) = Impulse4(X,t).*U(:,4);
   U(:,5) = Impulse4(X,t).*U(:,5);
   U(:,6) = Impulse5(X,t).*U(:,6);
   U(:,7) = Impulse5(X,t).*U(:,7);
   U(:,8) = Impulse6(X,t).*U(:,8);
   U(:,9) = Impulse6(X,t).*U(:,9);
end  
tcpu = cputime-tcpu;
fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
disp('Total Infectious Nymphs')
disp(sum(In(NN-1,:)))
%you'll need to change this if you change where the burn is being performed
disp('Infectious Nymphs in Burned Area')
disp(sum(In(NN-1,0 <= X(:,1) & X(:,1) <= 10 & 0 <= X(:,2) & X(:,2) <= 10)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PDE definition
function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    global piM miuM betaM piT sigmaT miuT betaT tauT gammaT miue K  
    Nm = u(:,8) + u(:,9);
    eps = 1e-10;
    lambdaT = betaT.*u(:,9)./(abs(Nm)+eps);
    lambdaM = betaM.*(u(:,3)+u(:,5)+u(:,7))./(abs(Nm)+eps);
    %diffusion functions
    DMX = Dm(x, t);
    DAX = Da(x, t);
    %pde model
    switch(ipde)
        case 1
        F = ut(:,1).*v(:) - piT.*(1-(u(:,1)/K)).*(u(:,6)+u(:,7)).*v(:) ...
            + (sigmaT + miue).*u(:,1).*v(:);
        case 2 %
        F = ut(:,2).*v(:) + DMX.*du(:,3).*dv(:,1)+ DMX.*du(:,4).*dv(:,2) - sigmaT.*u(:,1).*v(:) + lambdaT.*u(:,2).*v(:)...
            + (tauT+miuT).*u(:,2).*v(:);
        case 3 % 
        F = ut(:,3).*v(:) + DMX.*du(:,5).*dv(:,1)+ DMX.*du(:,6).*dv(:,2) - lambdaT.*u(:,2).*v(:) + (tauT+miuT).*u(:,3).*v(:);     
        case 4
        F = ut(:,4).*v(:) + DMX.*du(:,7).*dv(:,1)+ DMX.*du(:,8).*dv(:,2) - tauT.*u(:,2).*v(:) + lambdaT.*u(:,4).*v(:) ...
            + (gammaT+miuT).*u(:,4).*v(:);
        case 5 
        F = ut(:,5).*v(:) + DMX.*du(:,9).*dv(:,1)+ DMX.*du(:,10).*dv(:,2) - tauT.*u(:,3).*v(:) - lambdaT.*u(:,4).*v(:) ...
            + (gammaT+miuT).*u(:,5).*v(:);
        case 6 %
        F = ut(:,6).*v(:)+ DAX.*du(:,11).*dv(:,1)+ DAX.*du(:,12).*dv(:,2)  - gammaT.*u(:,4).*v(:) + lambdaT.*u(:,6).*v(:)...
            + miuT.*u(:,6).*v(:);
        case 7 %
        F = ut(:,7).*v(:)+ DAX.*du(:,13).*dv(:,1)+ DAX.*du(:,14).*dv(:,2)  - gammaT.*u(:,5).*v(:) - lambdaT.*u(:,6).*v(:)...
            + miuT.*u(:,7).*v(:);
        case 8 %
        F = ut(:,8).*v(:) + DMX.*du(:,15).*dv(:,1) + DMX.*du(:,16).*dv(:,2) - piM.*v(:) + lambdaM.*u(:,8).*v(:) + miuM.*u(:,8).*v(:); 
        case 9 % 
        F = ut(:,9).*v(:) + DMX.*du(:,17).*dv(:,1) + DMX.*du(:,18).*dv(:,2) - lambdaM.*u(:,8).*v(:) + miuM.*u(:,9).*v(:);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,ipde);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial Conditions
function u = u_initial(X)

   u = zeros(size(X,1),9);
%susceptible hosts present in new area:
   u(:,8) = 1000;
%    non uniform initial tick distribution
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,1) = 10000;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,2) = 0;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,3) = 100;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,4) = 0;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,5) = 500;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,6) = 0;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,7) = 100;
%    u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,8) = 1000;
   u(0 <= X(:,1) & X(:,1) <= 5 & 0<= X(:,2) & X(:,2) <=5,9) = 10;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diffusion functions
function D = Dm(X, t)

% both X and D are in dimesnional form

   global DM
   global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2

   D = zeros(size(X,1),1);
   D(:,:) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Da(X, t)

% both X and D are in dimesnional form

   global DA
   global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2

   D = zeros(size(X,1),1); 
   D(:,:) = 4;

%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse/prescribed fire functions
function Fe = Impulse1(X,t)
global nuE
global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2
% X is in dimesnional form, Fe is a dimensionless proportion
   Fe = ones(size(X,1),1);
   Fe(x1loc1 <= X(:,1) & X(:,1) <= x1loc2 & y1loc1 <= X(:,2) & X(:,2) <= y1loc2) = 1-nuE;
   Fe(x2loc1 <= X(:,1) & X(:,1) <= x2loc2 & y2loc1 <= X(:,2) & X(:,2) <= y2loc2) = 1-nuE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Flm = Impulse3(X,t)
global nuL
global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2
% X is in dimesnional form, Flm is a dimensionless proportion
   Flm = ones(size(X,1),1);
   Flm(x1loc1 <= X(:,1) & X(:,1) <= x1loc2 & y1loc1 <= X(:,2) & X(:,2) <= y1loc2) = 1-nuL;
   Flm(x2loc1 <= X(:,1) & X(:,1) <= x2loc2 & y2loc1 <= X(:,2) & X(:,2) <= y2loc2) = 1-nuL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fn = Impulse4(X,t)
global nuN
global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2
% X is in dimesnional form, Fn is a dimensionless proportion
   Fn = ones(size(X,1),1);  
   %Fn(0 <= X(:,1) & X(:,1) <= 0.5) = 1-nuN;
   Fn(x1loc1 <= X(:,1) & X(:,1) <= x1loc2 & y1loc1 <= X(:,2) & X(:,2) <= y1loc2) = 1-nuN;
   Fn(x2loc1 <= X(:,1) & X(:,1) <= x2loc2 & y2loc1 <= X(:,2) & X(:,2) <= y2loc2) = 1-nuN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fa = Impulse5(X,t)
global nuA
global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2
% X is in dimesnional form, Fa is a dimensionless proportion
   Fa = ones(size(X,1),1);  
   Fa(x1loc1 <= X(:,1) & X(:,1) <= x1loc2 & y1loc1 <= X(:,2) & X(:,2) <= y1loc2) = 1-nuA;
   Fa(x2loc1 <= X(:,1) & X(:,1) <= x2loc2 & y2loc1 <= X(:,2) & X(:,2) <= y2loc2) = 1-nuA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fm = Impulse6(X,t)
global nuM
global x1loc1 x1loc2 y1loc1 y1loc2 x2loc1 x2loc2 y2loc1 y2loc2
% X is in dimesnional form, Fm is a dimensionless proportion
   Fm = ones(size(X,1),1);  
   Fm(x1loc1 <= X(:,1) & X(:,1) <= x1loc2 & y1loc1 <= X(:,2) & X(:,2) <= y1loc2) = 1-nuM;
   Fm(x2loc1 <= X(:,1) & X(:,1) <= x2loc2 & y2loc1 <= X(:,2) & X(:,2) <= y2loc2) = 1-nuM;