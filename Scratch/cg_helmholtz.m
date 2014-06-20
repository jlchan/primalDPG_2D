% Poisson example for reference.
% square domain, forcing = 1
% BCs: left Dirichlet, right (nonconstant) inhomogeneous Neumann BC, and
% penalty enforcement of nonzero BCs on top and bottom

function cg_helmholtz

Globals2D

% Polynomial order used for approximation
N = 2;

global k
k = 1;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu'); % 8 elements
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu'); % 16 elements
% Initialize solver and construct grid and metric
StartUp2D;


PPB = numel(vmapB)/4 % 4 boundaries - points per boundary
PPW = 4;
k = PPB/PPW
% keyboard
% PPW = PPB/k

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);

% make CG operators
[R vmapBT] = getCGRestriction();

% forcing = 0
f = zeros(Np*K,1);
b = R*M*f;
B = R*BK*R';

% BC data for u
u0 = zeros(size(B,2),1);
t = 0*pi/4;
u0f = @(x,y) exp(1i*k*(cos(t)*x + sin(t)*y));
dudx0f = @(x,y) 1i*k*cos(t)*exp(1i*k*(cos(t)*x + sin(t)*y));
dudy0f = @(x,y) 1i*k*sin(t)*exp(1i*k*(cos(t)*x + sin(t)*y));

%%
% compute projections of exact solutions for bdata
Corder = 25;
[cubR,cubS,cubW, Ncub] = Cubature2D(Corder); Vcub = Vandermonde2D(N,cubR,cubS);
xcub = 0.5*(-(cubR+cubS)*VX(va)+(1+cubR)*VX(vb)+(1+cubS)*VX(vc));
ycub = 0.5*(-(cubR+cubS)*VY(va)+(1+cubR)*VY(vb)+(1+cubS)*VY(vc));

Interp = Vcub*invV; % interp to cubature points
Wcub = diag(cubW);
Minv = V*V';
u0 = Minv*Interp'*Wcub*u0f(xcub,ycub);
dudx0 = Minv*Interp'*Wcub*dudx0f(xcub,ycub);
dudy0 = Minv*Interp'*Wcub*dudy0f(xcub,ycub);

u0 = u0f(x(:),y(:));
dudx0 = dudx0f(x(:),y(:));
dudy0 = dudy0f(x(:),y(:));
dudn0 = nx(mapB).*dudx0(vmapB) + ny(mapB).*dudy0(vmapB);

%%

% penalty/robin BCs
% robin = y(vmapB) > 1 - NODETOL | x(vmapB) > 1-NODETOL; % top/right boundaries
robin = abs(nx(mapB)+1) < NODETOL | abs(ny(mapB)+1) < NODETOL;
% robin = zeros(size(vmapB));
[Mb Eb] = getBoundaryMatrix(robin);
B = B - 1i*k*R*Eb'*Mb*Eb*R';
b = b + R*Eb'*Mb*(-1i*k*u0(vmapB) + dudn0);

% nonhomogeneous neumann BCs
% onlyLeft = (x(vmapB) < -1 + NODETOL) & (y(vmapB) < 1-NODETOL);
% onlyBot = (y(vmapB) < -1 + NODETOL) & (x(vmapB) < 1-NODETOL);
% neum = onlyLeft | onlyBot; % left/bottom boundary
neum = abs(nx(mapB)-1) < NODETOL | abs(ny(mapB)-1) < NODETOL;
% neum = ones(size(vmapB))*0;
[Mb Eb] = getBoundaryMatrix(neum(:));
b = b + R*Eb'*Mb*dudn0;

% solve and prolong solution u to local storage
U = (B\b);
u = R'*U;

% Nplot = 25; [xu,yu] = EquiNodes2D(Nplot);
% Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot);
Nplot = 25; [xu,yu] = Nodes2D(Nplot);
[ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));

% figure
Iu0 = u0f(xu(:),yu(:));

figure
Iuh = Iu*reshape(u,Np,K);Iuh = Iuh(:);
ax = [-1 1 -1 1 -1 1];
subplot(3,1,1);color_line3(xu,yu,Iuh,Iuh,'.');axis(ax);view(2);%view(0,0)
subplot(3,1,2);color_line3(xu,yu,Iu0,Iu0,'.');axis(ax);view(2);%view(0,0)
subplot(3,1,3);color_line3(xu,yu,Iuh-Iu0,Iuh-Iu0,'.');%axis([-1 1 -1 1 -.1 .1])
colorbar;view(2);%view(0,0);

Vu = Vandermonde2D(Nplot,ru,su);
Mu = inv(Vu*Vu');
Nplotp = (Nplot+1)*(Nplot+2)/2;
e = reshape(Iuh-Iu0,Nplotp,K);
J = J(1,:); % hack for non curved elems
J = repmat(J,Nplotp,1);
tmp = J.*(Mu*e);
L2err = sqrt(abs(e(:)'*tmp(:)));
title(['L2 error = ' num2str(L2err) ' at k = ' num2str(k) ' and ' num2str(PPW) ' ppw.'])

% keyboard

function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
% Test = M + Ks;
% Trial = Ks;

global k;
Trial = -k^2*M + Ks;
% Test = k^2*M + Ks;
Test = k^2*M + Ks;

% blkiM = kron(speye(K),inv(MassMatrix));
% invM = spdiag(1./J(:))*blkiM; % J = h^2
% Test = Trial*blkiM*Trial' + 1e-6*M;%(1/k^2)*M;

