function [A b] = UWDPG_poisson

Globals2D
FaceGlobals2D;

% Polynomial order used for approximation
Ntrial = 2;
Ntest = Ntrial + 2;
Nf = Ntrial;
Nt = Nf;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
[Nv, VX, VY, K, EToV] = QuadMesh2D(4);
% [Nv, VX, VY, K, EToV] = MakeQuads2D(4);
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('backdrop1.neu');

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);nyfb = nyf(fmapB);

% get block operators
[M, Dx, Dy] = getBlockOps();
M2 = blkdiag(M,M); O = zeros(size(M));
Div = [Dx Dy];
Grad = [Dx;Dy];
AK = blkdiag(M2 + Div'*M*Div,M + Grad'*M2*Grad);
BK = [M2 Div'*M;Grad'*M2 O];
f = y(:).^0;

Bhat1 = getMortarConstraint();
Bhat2 = getMortarConstraintDiv();
Bhat = blkdiag(Bhat2',Bhat1');

NpTest = Np;

% switches to Ntrial globals
[~, Rr] = pRestrictCG(N,Ntrial); % restrict test to trial space
B = BK*blkdiag(Rr,Rr,Rr);   % form rectangular bilinear form matrix
[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,2); % num mortar nodes
% nTrial = nU + nM;

Bh = [B Bhat];
Tblk = cell(K,1);
p = [];
for e = 1:K % independently invert
    inds1 = (e-1)*NpTest + (1:NpTest); % tau1
    inds2 = NpTest*K + (e-1)*NpTest + (1:NpTest); % tau2
    inds3 = 2*NpTest*K + (e-1)*NpTest + (1:NpTest); % v
    inds = [inds1 inds2 inds3];    
    Tblk{e} = AK(inds,inds)\Bh(inds,:);
    disp(['on element ' num2str(e)])
    p = [p inds]; % permutation inds
end
T = cell2mat(Tblk);
T(p,:) = T;
A = T'*Bh;

% forcing
b = T'*blkdiag(M2,M)*[zeros(2*NpTest*K,1);f];

% BCs on u
u0 = zeros(nU,1);

% BCs on flux
uh0 = zeros(nM,1);
U0 = [u0;uh0];

% BCs on U: ordered first
b = b - A*U0;

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
bci = nU + tmapB; % skip over u dofs
b(bci) = 0;
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

U = A\b;
u = U(2*Np*K+1:3*Np*K);

Nplot = 25;
plotSol(u,Nplot);
% plotFlux(f)
title('DPG with fluxes and traces')


