function mortarTestNorm

Globals2D;
FaceGlobals2D

N = 7; % when N = even, Nf = N-1, fails?
Nf = 6; % = N trial
Nt = 7; % = 
%     Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
%     Nv = 3;
%     VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
%     EToV = [3 1 2];
%     K = 1;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');

% Initialize solver and construct grid and metric
StartUp2D; FaceStartUp2D

[M, Dx, Dy] = getBlockOps();
Div = [Dx Dy];
Grad = [Dx;Dy];
M2 = blkdiag(M,M);
I2 = speye(size(M2));
I = speye(size(M));
O = zeros(size(M));

% % Poisson
% Adj_h = [I2 Grad;
%          Div O];

% % convection-diffusion
% ep = .1;
% Adj_h = [(1/ep)*I2 Grad;
%          Div -Dx];

% Helmholtz
k = 25;
Adj_h = [1i*k*I2 Grad;
         Div 1i*k*I];

M3 = blkdiag(M2,M);
RV = Adj_h'*M3*Adj_h + 1e-4*blkdiag(M,M,M); % regularize on v

f = x(:).^0;
btau = M2*[f;f];
bv = M*f*0;
% bv = bv*0;
% bv(end-10) = 1;
b = Adj_h'*[btau;bv];

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
B = [Bt zeros(mtau,nv);
     zeros(mv,ntau) Bf];
 
nU = size(B,2); % num CG nodes
nM = size(B,1); % num mortar nodes

Am = [RV B'
    B zeros(nM)];
bm = [b;zeros(nM,1)];
um = Am\bm;

U = um(1:nU);

I1 = 1:Np*K;
I2 = Np*K + (1:Np*K);
I3 = 2*Np*K + (1:Np*K);

tau1 = U(I1);
tau2 = U(I2);
v = U(I3);

taunorm = sqrt(tau1.^2 + tau2.^2);
plotSol(taunorm,25);title('norm of tau')
plotSol(v,25);title('v')
