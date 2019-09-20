function top3de(nelx,nely,nelz,volfrac,penal,rmin,ft,extdirection)
%Example inputs: nelx=60;nely=20;nelz=20;volfrac=0.3;penal=3;rmin=3.5;extdirection=3;ft=3;
%% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
%% USER-DEFINED MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-6;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio
%% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
%% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
%% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8s(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
%% PREPARE FILTER
Rmin = rmin/2/sqrt(3);
c = Rmin^2*ones(1,3); 
c(extdirection) = Rmin^2*2500;%Filter Radius of the extrusion direction
KEF = lk_H8f(c);
edofVecF = nodeids(:)+1;
edofMatF = repmat(edofVecF,1,8)+ ...
    repmat([0 nely + [1 0] -1 ...
    (nely+1)*(nelx+1)+[0 nely + [1 0] -1]],nele,1);
iKF = reshape(kron(edofMatF,ones(8,1))',8*8*nele,1);
jKF = reshape(kron(edofMatF,ones(1,8))',8*8*nele,1);
sKF = reshape(KEF(:)*ones(1,nele),8*8*nele,1);
KF = sparse(iKF,jKF,sKF);
LF = chol(KF,'lower');
iTF = reshape(edofMatF,8*nele,1);
jTF = reshape(repmat([1:nele],8,1)',8*nele,1);
sTF = repmat(1/8,8*nele,1);
TF = sparse(iTF,jTF,sTF);
%% MMA SETTINGS
upp=1;
low=0;
m=1;
n=nele;
ccc=1000;
a0 = 1;
a = zeros(m,1);
d = ones(m,1);
%% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
xold1 = x;
xold2 = x;
beta = 2;
xPhys = x;
if ft == 3
   xTilde = x;
   eta=0.5;
   xPhys=(tanh(beta.*eta)+tanh(beta.*(xTilde-eta)))./(tanh(beta.*eta)+tanh(beta.*(1-eta)));
end
loop = 0; 
loopbeta = 0;
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    loopbeta = loopbeta+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx,nelz);
    %% FILTERING AND MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = (TF'*(LF'\(LF\(TF*(dc(:).*xPhys(:))))))./max(1e-3,xPhys(:));
    elseif ft == 2
        dc(:) = TF'*(LF'\(LF\(TF*dc(:))));
        dv(:) = TF'*(LF'\(LF\(TF*dv(:))));
    elseif ft == 3
        dx = -(beta.*(tanh(beta.*(eta - xTilde(:))).^2 - 1))./(tanh(beta.*eta) - tanh(beta.*(eta - 1)));
        dc(:) = TF'*(LF'\(LF\(TF*(dc(:).*dx(:)))));
        dv(:) = TF'*(LF'\(LF\(TF*(dv(:).*dx(:)))));
    end
    %% MMA UPDATE OF DESIGN VARIABLES
    [xnew,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,x(:),0,1,xold1(:),xold2(:), ...
        c,dc(:),sum(xPhys(:))-volfrac*nele,dv(:)',low,upp,a0,a,ccc,d);
    change = max(abs(xnew(:)-x(:)));
    xold1 = x;
    xold2 = xold1;
    x = xnew;
    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if ft == 3 && beta < 64 && (loopbeta >= 30 || change <= 0.01)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
    %% UPDATE PHYSICAL DENSITIES
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (TF'*(LF'\(LF\(TF*xnew(:)))));
    elseif ft == 3
        xTilde(:) = (TF'*(LF'\(LF\(TF*xnew(:)))));
        xPhys=(tanh(beta.*eta)+tanh(beta.*(xTilde-eta)))./(tanh(beta.*eta)+tanh(beta.*(1-eta)));
    end
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    %% PLOT DENSITIES
    if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
end
clf; display_3D(xPhys);
end

%% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8s(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end

%% === GENERATE ELEMENT HELMHOLTZ-OPERATOR MATRIX ===
function [HE] = lk_H8f(c)
H1= [   1/9,  -1/9, -1/18,  1/18,  1/18, -1/18, -1/36,  1/36;
       -1/9,   1/9,  1/18, -1/18, -1/18,  1/18,  1/36, -1/36;
      -1/18,  1/18,   1/9,  -1/9, -1/36,  1/36,  1/18, -1/18;
       1/18, -1/18,  -1/9,   1/9,  1/36, -1/36, -1/18,  1/18;
       1/18, -1/18, -1/36,  1/36,   1/9,  -1/9, -1/18,  1/18;
      -1/18,  1/18,  1/36, -1/36,  -1/9,   1/9,  1/18, -1/18;
      -1/36,  1/36,  1/18, -1/18, -1/18,  1/18,   1/9,  -1/9;
       1/36, -1/36, -1/18,  1/18,  1/18, -1/18,  -1/9,   1/9];
H2= [   1/9,  1/18, -1/18,  -1/9,  1/18,  1/36, -1/36, -1/18;
       1/18,   1/9,  -1/9, -1/18,  1/36,  1/18, -1/18, -1/36;
      -1/18,  -1/9,   1/9,  1/18, -1/36, -1/18,  1/18,  1/36;
       -1/9, -1/18,  1/18,   1/9, -1/18, -1/36,  1/36,  1/18;
       1/18,  1/36, -1/36, -1/18,   1/9,  1/18, -1/18,  -1/9;
       1/36,  1/18, -1/18, -1/36,  1/18,   1/9,  -1/9, -1/18;
      -1/36, -1/18,  1/18,  1/36, -1/18,  -1/9,   1/9,  1/18;
      -1/18, -1/36,  1/36,  1/18,  -1/9, -1/18,  1/18,   1/9];
H3= [   1/9,  1/18,  1/36,  1/18,  -1/9, -1/18, -1/36, -1/18;
       1/18,   1/9,  1/18,  1/36, -1/18,  -1/9, -1/18, -1/36;
       1/36,  1/18,   1/9,  1/18, -1/36, -1/18,  -1/9, -1/18;
       1/18,  1/36,  1/18,   1/9, -1/18, -1/36, -1/18,  -1/9;
       -1/9, -1/18, -1/36, -1/18,   1/9,  1/18,  1/36,  1/18;
      -1/18,  -1/9, -1/18, -1/36,  1/18,   1/9,  1/18,  1/36;
      -1/36, -1/18,  -1/9, -1/18,  1/36,  1/18,   1/9,  1/18;
      -1/18, -1/36, -1/18,  -1/9,  1/18,  1/36,  1/18,   1/9];
H4= [  1/27,  1/54, 1/108,  1/54,  1/54, 1/108, 1/216, 1/108;
       1/54,  1/27,  1/54, 1/108, 1/108,  1/54, 1/108, 1/216;
      1/108,  1/54,  1/27,  1/54, 1/216, 1/108,  1/54, 1/108;
       1/54, 1/108,  1/54,  1/27, 1/108, 1/216, 1/108,  1/54;
       1/54, 1/108, 1/216, 1/108,  1/27,  1/54, 1/108,  1/54;
      1/108,  1/54, 1/108, 1/216,  1/54,  1/27,  1/54, 1/108;
      1/216, 1/108,  1/54, 1/108, 1/108,  1/54,  1/27,  1/54;
      1/108, 1/216, 1/108,  1/54,  1/54, 1/108,  1/54,  1/27];
HE=c(1)*H1+c(2)*H2+c(3)*H3+H4;
end

%% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Zhou Yan, Department of Engineering      %
% Mechanics, State Key Laboratory of Structural Analysis for Industrial    %
% Equipment                                                                %
% Dalian University of Technology                                          %
% Dalian, China                                                            %
%                                                                          %
% The code is intended for replication of results in the papaer            %
% "Novel implementation of extrusion constraint in topology optimization   %
%  by Helmholtz-type anisotropic filter".                                  %
%                                                                          %
% This code is modified from the following codes include                   %
% TOP88, TOP82 in "Andreassen, Erik, et al. Efficient topology             %
% optimization in MATLAB using 88 lines of code.                           %
% Structural and Multidisciplinary Optimization 43.1 (2011): 1-16."        %
% and                                                                      %
% TOP3D in "Liu, Kai, and Andrés Tovar. An efficient 3D topology           %
% optimization code written in Matlab.                                     %
% Structural and Multidisciplinary Optimization 50.6 (2014): 1175-1196."   %                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%