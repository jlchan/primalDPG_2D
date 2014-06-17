% Poisson example for reference. 
% square domain, forcing = 1
% BCs: left Dirichlet, right (nonconstant) inhomogeneous Neumann BC, and
% penalty enforcement of nonzero BCs on top and bottom 

function mixedDPG_poisson

Globals2D

% Polynomial order used for approximation
Ntrial = 4;
Ntest = Ntrial+0;

N = Ntest;
global k
k = 25;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);

% forcing = 1
f = 0*ones(Np*K,1);

% make CG operators
[R vmapBT] = getCGRestriction();
[Rr vmapBTr xr yr] = pRestrictCG(Ntrial); % restrict test to trial space 
B = R*BK*Rr';
RV = R*AK*R';
b = R*M*f;

% BC data for u
u0 = zeros(size(B,2),1);
t = pi/4;
u0f = @(x,y) exp(1i*k*(cos(t)*x + sin(t)*y));
u0 = u0f(x(:),y(:));
dudx0 = 1i*k*cos(t)*exp(1i*k*(cos(t)*x(:) + sin(t)*y(:)));
dudy0 = 1i*k*sin(t)*exp(1i*k*(cos(t)*x(:) + sin(t)*y(:)));
dudn0 = nx(mapB).*dudx0(vmapB) + ny(mapB).*dudy0(vmapB);

% penalty/robin BCs 
bmask = y(vmapB) > 1 - NODETOL | x(vmapB) > 1-NODETOL; % top/right boundaries
[Mb Eb] = getBoundaryMatrix(bmask);
B = B + 1i*k*R*Eb'*Mb*Eb*Rr';
b = b + R*Eb'*Mb*(1i*k*u0(vmapB) + dudn0);

% t = 0;
% uex = exp(i*k*(cos(t)*x(:)+sin(t)*y(:)));
% g = i*k*uex(vmapB);
% b = b + sqrt(-1)*k*R*Eb'*Mb*g;
% B = B + sqrt(-1)*k*R*Eb'*Mb*Eb*Rr';

% nonhomogeneous neumann BCs
onlyLeft = (x(vmapB) < -1 + NODETOL) & (y(vmapB) < 1-NODETOL);
onlyBot = (y(vmapB) < -1 + NODETOL) & (x(vmapB) < 1-NODETOL);
bmask = onlyLeft | onlyBot; % left/bottom boundary
[Mb Eb] = getBoundaryMatrix(bmask(:)); 
b = b + R*Eb'*Mb*dudn0;

% make saddle point system
A = [RV B;B' zeros(size(B,2))];
b = [b; zeros(size(B,2),1)]; 

% solve and prolong solution u to local storage
U = (A\b);
u = Rr'*U(size(B,1)+1:end);

e = U(1:size(B,1));
% err = e'*RV*e;
e = R'*e;

Nplot = 25; [xu,yu] = EquiNodes2D(Nplot); 
% Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot); 
[ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));

% figure
Iu0 = u0f(xu(:),yu(:));
% color_line3(xu,yu,Iu0,Iu0,'.');

figure 
Iuh = Iu*reshape(u,Np,K);Iuh = Iuh(:);
Ie = Iu*reshape(e,Np,K);Ie = Ie(:);
title('Mixed form of DPG')
ax = [-1 1 -1 1 -1 1];
subplot(2,1,1);color_line3(xu,yu,Iuh,Iuh,'.');axis(ax);view(2);%view(0,0)
subplot(2,1,2);color_line3(xu,yu,Iu0,Iu0,'.');axis(ax);view(2);%view(0,0)
figure
subplot(2,1,1);color_line3(xu,yu,Ie,Ie,'.');axis([-1 1 -1 1 -.1 .1])
colorbar;view(2);%view(0,0);
subplot(2,1,2);color_line3(xu,yu,Iuh-Iu0,Iuh-Iu0,'.');axis([-1 1 -1 1 -.1 .1])
colorbar;view(2);%view(0,0);


% keyboard

function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
% Test = M + Ks;
% Trial = Ks;

global k;
Trial = -k^2*M + Ks;
Test = k^2*M + Ks;
% blkiM = kron(speye(K),inv(MassMatrix));
% invM = spdiag(1./J(:))*blkiM; % J = h^2
% Test = Trial*blkiM*Trial' + 1e-6*M;%(1/k^2)*M;

function [M, Dx, Dy] = getBlockOps()

Globals2D

blkDr = kron(speye(K),Dr);
blkDs = kron(speye(K),Ds);
blkM = kron(speye(K),MassMatrix);

M = spdiag(J(:))*blkM; % J = h^2
Dx = spdiag(rx(:))*blkDr + spdiag(sx(:))*blkDs;
Dy = spdiag(ry(:))*blkDr + spdiag(sy(:))*blkDs;


