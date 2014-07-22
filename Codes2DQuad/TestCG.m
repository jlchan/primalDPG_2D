clear
% Driver script for solving the 2D vacuum Maxwell's equations on TM form
% Globals2D;
Globals2D;

% Polynomial order used for approximation
N = 3;

[Nv, VX, VY, K, EToV] = QuadMesh2D(5);

StartUp2D;

[M, Dx, Dy] = getBlockOps();
[R vmapBT] = getCGRestriction();
Grad = [Dx;Dy];
A = M + Grad'*blkdiag(M,M)*Grad;

A = R*A*R';
f = ones(Np*K,1);
b = R*M*f;

n = size(A,1);
b(vmapBT) = 0;
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));
u = A\b;
u = R'*u;

color_line3(x,y,u,u,'.')
plotSol(u,25)