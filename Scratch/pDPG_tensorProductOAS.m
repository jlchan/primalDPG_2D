function pDPG_tensorProductOAS

%% make DPG system
Globals2D
FaceGlobals2D

% Polynomial order used for approximation
Ntrial = 2;
Ntest = 4;
Nf = 2;
% Nt = 0; % dummy

N = Ntest;
[Nv, VX, VY, K, EToV] = QuadMesh2D(2);

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D 

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = ones(Np*K,1);
b = M*f;

%[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xrB yrB xr yr] = pRestrictCG(Ntest,Ntrial); % restrict test to trial space
Rr = Rp*Irp';
Bhat = getMortarConstraint();
xfB = xf(fmapB); yfB = yf(fmapB); nxf = nxf(fmapB);nyf = nyf(fmapB);

B = BK*Rr';   % form rectangular bilinear form matrix

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
tic
for i = 1:K % independently invert
    inds = (i-1)*Np + (1:Np);
    Tblk{i} = AK(inds,inds)\Bh(inds,:);
    disp(['on element ' num2str(i)])    
end
disp(['time for test function computation = ', num2str(toc)])
T = cell2mat(Tblk);

% assemble forcing + matrix
A = T'*Bh;
b = T'*b;

% BCs on u, f
u0 = zeros(size(B,2),1);
uh0 = zeros(nM,1); 

U0 = [u0;uh0];
b = b - A*U0; % get lift

% BCs on U: ordered first
b(vmapBTr) = U0(vmapBTr); 
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0; A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

U = A\b;

u = Rr'*U(1:nU);
uhat = U(nU+(1:nM));

figure
plotFlux(uhat)
plotSol(u,25)
% plotSol(Dx*u,25)
% plotSol(Dy*u,25)

title('DPG with fluxes and traces')

%% build OAS dof lists
[Mhandle Ak] = buildOAS_primalDPG(Rp,A,Ntrial,2); % 2 = face patch, 3 = elem patch only

NpU = (Ntrial+1)^2;
NpF = Nfrp*4;
iU = 1:NpU;
iF = NpU+1:NpU+NpF; 

keyboard
function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = M + Ks;
Trial = Ks;

