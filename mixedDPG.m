function mixedDPG

Globals2D

% Polynomial order used for approximation
Ntrial = 3;
Ntest = Ntrial + 2;

N = Ntest;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% % [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('backdrop1.neu');

% Initialize solver and construct grid and metric
StartUp2D;

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

[R vmapBT] = getCGRestriction();
[Rr vmapBTr xr yr] = pRestrictCG(Ntrial); % restrict test to trial space 

% make CG operators
B = R*BK*Rr';
RV = R*AK*R';
b = R*M*f;

% BC data for u
u0 = zeros(size(B,2),1);
% u0(vmapBTr) = (xr < -1+1e-7).*sqrt(1-yr.^2); 
% to fix - add integral over boundaries for inhomog Neumann conditions
bnf = nx(mapB)*b1 + ny(mapB)*b2; % beta_n, determines inflow vs outflow
bmask = (bnf < NODETOL); % inflow = beta_n < 0
f0 = bnf.*(Fy(mapB)<0).*(1+Fy(mapB));
[Mb Eb] = getBoundaryMatrix(bmask.^0);
b = b - R*Eb'*Mb*f0;  % BC data on flux = bn*u - eps*du/dn

% BC data for e is generally zero.  
e0 = zeros(size(B,1),1);
b = [b; zeros(size(B,2),1)]; 
U0 = [e0;u0];

% B = B + 1e7*R*Eb'*Mb*Eb*Rr'; % this adds a penalty term

% make saddle point system
A = [RV B;B' zeros(size(B,2))];
% A = [RV B;B' 1e2*Rr*Mb*Rr']; %penalty BCs as functional

% applying lift data
b = b - A*U0;

% BCs on U - skip over e dofs
outflow = xr>1-NODETOL;
vmapBTr(~outflow) = [];
vmapBTU = vmapBTr + size(B,1);
b(vmapBTU) = U0(vmapBTU);
A(vmapBTU,:) = 0; A(:,vmapBTU) = 0;
A(vmapBTU,vmapBTU) = speye(length(vmapBTU));

% BCs on V
beta_n = b1*nx + b2*ny; beta_nb = beta_n(mapB);
inflow = beta_nb < -NODETOL;
vmapBT(inflow)=[]; % remove BCs on inflow for convection-diffusion!
% vmapBT= []; %remove all BCs for testing
b(vmapBT) = U0(vmapBT);
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));

% solve and prolong solution u to local storage
U = (A\b);
u = Rr'*U(size(B,1)+1:end);
% e = U(1:size(B,1));
% err = e'*RV*e;
% e = R'*e;

Nplot = 25;
[xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
title('Mixed form of DPG')

% keyboard

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
Test = M + Ks;
Trial = Ks;

% CD
Test = M + ep*Ks + Kb;
Trial = S + ep*Ks;

% Helmholtz
% k = 100;
% Test = k^2*M + Ks;
% Trial = -k^2*M + Ks;

function [M, Dx, Dy] = getBlockOps()

Globals2D

blkDr = kron(speye(K),Dr);
blkDs = kron(speye(K),Ds);
blkM = kron(speye(K),MassMatrix);

M = spdiag(J(:))*blkM; % J = h^2
Dx = spdiag(rx(:))*blkDr + spdiag(sx(:))*blkDs;
Dy = spdiag(ry(:))*blkDr + spdiag(sy(:))*blkDs;


