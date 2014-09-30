function primalDPG_confusion

Globals2D
FaceGlobals2D;

% Polynomial order used for approximation
Ntrial = 6;
Ntest = Ntrial + 2;
Nf = Ntrial;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
[Nv, VX, VY, K, EToV] = QuadMesh2D(3);
% [Nv, VX, VY, K, EToV] = MakeQuads2D(4);
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('backdrop1.neu');

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

global b1
global b2
global ep
b1 = 1; b2 = 0;ep = 1e-6;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = 0*ones(Np*K,1);
% f = y(:)<=0;
% f = sin(pi*x(:)).*sin(pi*y(:));

Bhat = getMortarConstraint();
xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);nyfb = nyf(fmapB);

NpTest = Np;
[Rp Irp vmapBT] = pRestrictCG(N,Ntrial); % restrict test to trial space
Rr = Rp*Irp';
B = BK*Rr';   % form rectangular bilinear form matrix
[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
if 1
    tic
    for i = 1:K % independently invert
        inds = (i-1)*NpTest + (1:NpTest);
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
        disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
    disp(['time for test function computation = ', num2str(toc)])
else
    disp('parfor implementation...')
    t = 0;
    tic
    AKi = cell(K,1); Bi = cell(K,1);
    for i = 1:K
        inds = (i-1)*Np + (1:Np);
        AKi{i} = AK(inds,inds);
        Bi{i} = Bh(inds,:);
    end
    t=t+toc
    tic
    parfor i = 1:K
        Tblk{i} = AKi{i}\Bi{i};
    end
    t = t+toc
end
T = cell2mat(Tblk);
A = T'*Bh;

% Grad = [Dx;Dy];
% Mxy = M;
% reg = Grad'*blkdiag(Mxy,Mxy)*Grad;
% h = spdiag(J(:));
% reg = M;
% A(1:nU,1:nU) = A(1:nU,1:nU) + Rr*h*reg*Rr';

% forcing
b = T'*M*f;

% BCs on u
u0 = zeros(size(B,2),1);

% BCs on flux
uh0 = zeros(nM,1);
bnf = nxfb*b1 + nyfb*b2; % beta_n, determines inflow vs outflow
bmaskf = (bnf < NODETOL); % inflow = beta_n < 0
uh0(fmapB) = bnf.*(yfb<0).*(1+yfb);  % BC data on flux = bn*u - eps*du/dn
% uh0(fmapB) = bnf;  % BC data on flux = bn*u - eps*du/dn

U0 = [u0;uh0];

% remove BCs on u on inflow for stability
vmapBT(x(vmapB) < -1+NODETOL) = [];
% wall = (abs(yr+1)<NODETOL) & (xr > -NODETOL);
% vmapBTr(~wall) = [];
% vmapBT = []; % removes all Dirichlet BCs for testing....

% BCs on U: ordered first
b = b - A*U0;
b(vmapBT) = U0(vmapBT);
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
fmapB(bmaskf) = []; % do 0 Neumann outflow BCs on test fxns
% fmapB = [];
% fmapB = fmapB(bmaskf);

bci = nU + fmapB; % skip over u dofs
b(bci) = uh0(fmapB);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

U = A\b;
u = Rp'*U(1:nU);
f = U(nU+(1:nM));
%     color_line3(x,y,u,u,'.');
%     return

Nplot = 25;
plotSol(u,Nplot);
plotFlux(f)
title('DPG with fluxes and traces')




function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
global b1
global b2
global ep
Ks = Dx'*M*Dx + Dy'*M*Dy;

S = -(b1*Dx+b2*Dy)'*M;
% S = M*(b1*Dx+b2*Dy);
Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Poisson
% Test = M + Ks;
% Trial = Ks;

% CD
Test = M + ep*Ks + Kb;
Trial = S + ep*Ks;

% Helmholtz
% k = 100;
% Test = k^2*M + Ks;
% Trial = -k^2*M + Ks;

