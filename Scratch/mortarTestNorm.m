function mortarTestNorm

Globals2D;
FaceGlobals2D

N = 6; % when N = even, Nf = N-1, fails?
Nf = 2; % = N trial
Nt = 3; % = 
%     Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
%     Nv = 3;
%     VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
%     EToV = [3 1 2];
%     K = 1;
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu'); 
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');

[Nv, VX, VY, K, EToV] = QuadMesh2D(4);
% [Nv, VX, VY, K, EToV] = MakeQuads2D(8);

% Initialize solver and construct grid and metric
StartUp2D; FaceStartUp2D

[M, Dx, Dy] = getBlockOps();
Div = [Dx Dy];
Grad = [Dx;Dy];
M2 = blkdiag(M,M);
I2 = speye(size(M2));
I = speye(size(M));
O = sparse(size(M,1),size(M,2));

% % Poisson
% Adj_h = [I2 Grad;
%          Div O];

% % convection-diffusion
% ep = .01;
% Adj_h = [(1/ep)*I2 Grad;
%          Div -Dx];

% Helmholtz
% px = sqrt(K)*(N+1);
% k = 10;
% Adj_h = [1i*k*I2 Grad;
%          Div 1i*k*I];

M3 = blkdiag(M2,M);
RV = Adj_h'*M3*Adj_h + 0*blkdiag(M,M,M); % regularize on v

% time-dependent heat eqn
ep = .01;
Adj_h = [(1/ep)*I2 blkdiag(Dx);
         Div O];


f = x(:).^0;
btau = M2*[f;f];
bv = M*f;
% bv = bv*0;
% bv(end-10) = 1;
b = Adj_h'*M3*[btau;bv];

Bt = getMortarConstraintDiv();
% fmapBd = fmapB;xfd = xf;yfd = yf; 
xtb = xt(tmapB);ytb = yt(tmapB); nxtb = nxt(tmapB);
Bt(tmapB((abs(1-xtb)<NODETOL | abs(1-ytb.^2) < NODETOL) & nxtb > -NODETOL),:) = []; %remove constraints for fluxes
% Bt(tmapB,:) = [];

Bf = getMortarConstraint();
xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);
Bf(fmapB(abs(1+xfb)<NODETOL & abs(nxfb+1)<NODETOL),:) = []; %remove constraints for fluxes

[mtau ntau] = size(Bt);
[mv nv] = size(Bf);
B = [Bt sparse(mtau,nv);
     sparse(mv,ntau) Bf];
% B = B*Adj_h; % change trace space
nU = size(B,2); % num CG nodes
nM = size(B,1); % num mortar nodes

Am = [RV B'
    B sparse(nM,nM)];
bm = [b;zeros(nM,1)];
um = Am\bm;
U = um(1:nU);

% S = B*(RV\B');
% bS = B*(RV\b);
% f = S\bS;
% U = RV\(b - B'*f);

I1 = 1:Np*K;
I2 = Np*K + (1:Np*K);
I3 = 2*Np*K + (1:Np*K);

tau1 = U(I1);
tau2 = U(I2);
v = U(I3);

taunorm = sqrt(tau1.^2 + tau2.^2);
plotSol(taunorm,25);%title(['norm of tau, px = ', num2str(px)])
plotSol(v,25);%title(['v, px = ', num2str(px)])
if (Np*K<500)
    figure;semilogy(svd(full(Adj_h)))
end
% keyboard
