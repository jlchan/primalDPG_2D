% Poisson example for reference. 
% square domain, forcing = 1
% BCs: left Dirichlet, right (nonconstant) inhomogeneous Neumann BC, and
% penalty enforcement of nonzero BCs on top and bottom 

function mixedDPG_poisson

Globals2D

% Polynomial order used for approximation
Ntrial = 2;
Ntest = Ntrial+2;

N = Ntest;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);

% forcing = 1
f = ones(Np*K,1);

% make CG operators
[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xr yr] = pRestrictCG(N,Ntrial); % restrict test to trial space 
Rr = Rp*Irp';
B = R*BK*Rr';
RV = R*AK*R';
b = R*M*f;

% BC data for u
u0 = zeros(size(B,2),1);
left = xr < -1+NODETOL;
right = xr > 1-NODETOL;
u0(vmapBTr) = left.*sqrt(1-yr.^2); 

% penalty/robin BCs 
bmask = abs(y(vmapB)) > 1 - NODETOL; % top/bottom boundaries
[Mb Eb] = getBoundaryMatrix(bmask);
u0tb = 1+x(vmapB);
B = B + 1e6*R*Eb'*Mb*Eb*Rr'; % this adds a penalty term on u (or Robin condition) 
b = b + 1e6*R*Eb'*Mb*u0tb;

% nonhomogeneous neumann BCs
bmask = x(vmapB) > 1 - NODETOL; % right boundary
[Mb Eb] = getBoundaryMatrix(bmask(:)); 
dudn0 = nx(mapB).*((y(vmapB)<=0) - (y(vmapB)>0));
b = b + R*Eb'*Mb*dudn0;

% BC data for e is generally zero.  
e0 = zeros(size(B,1),1);
b = [b; zeros(size(B,2),1)]; 
U0 = [e0;u0];

% make saddle point system
A = [RV B;B' zeros(size(B,2))];

% applying lift data
b = b - A*U0;

% BCs on U - skip over e dofs
vmapBTr(~left) = []; %neumann on right outflow - remove top/bottom BCs
% vmapBTr(right) = []; %neumann on right outflow - remove top/bottom BCs
vmapBTU = vmapBTr + size(B,1);
b(vmapBTU) = U0(vmapBTU);
A(vmapBTU,:) = 0; A(:,vmapBTU) = 0;
A(vmapBTU,vmapBTU) = speye(length(vmapBTU));

% BCs on V
left = x(vmapB)<-1+NODETOL; 
right = x(vmapB) > 1-NODETOL;
vmapBT(~left) = []; % remove BCs for dirichlet conditions only
b(vmapBT) = U0(vmapBT);
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));

% solve and prolong solution u to local storage
U = (A\b);
u = Rr'*U(size(B,1)+1:end);
% e = U(1:size(B,1));
% err = e'*RV*e;
% e = R'*e;

Nplot = 25; [xu,yu] = EquiNodes2D(Nplot); 
Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot); 
[ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
title('Mixed form of DPG')

% keyboard

function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = M + Ks;
Trial = Ks;

