% Poisson example for reference.
% square domain, forcing = 1
% BCs: left Dirichlet, right (nonconstant) inhomogeneous Neumann BC, and
% penalty enforcement of nonzero BCs on top and bottom

function cg_helmholtz

Globals2D

% Polynomial order used for approximation
N = 10;

global k
k = 1;

[Nv, VX, VY, K, EToV] = QuadMesh2D(1);
StartUp2D;

k = 5;

% get block operators
M = spdiag(J(:))*kron(speye(K),inv(V*V')); % J = h^2
Dx = spdiag(rx(:))*kron(speye(K),Dr) + spdiag(sx(:))*kron(speye(K),Ds);
Dy = spdiag(ry(:))*kron(speye(K),Dr) + spdiag(sy(:))*kron(speye(K),Ds);

Ks = Dx'*M*Dx + Dy'*M*Dy;
A = -k^2*M + Ks;

% forcing = 0
f = zeros(Np,1);
b = M*f;

% BC data for u
t = 0*pi/4;
[u Dxu Dyu] = exact_sol_PW(x(vmapB),y(vmapB),k,0);
dudn = nx(mapB).*Dxu + ny(mapB).*Dyu;
g = dudn - 1i*k*u;
    
% robin BCs
[Mb Eb] = getBoundaryMatrix();
A = A - 1i*k*Eb'*Mb*Eb; 
b = b + Eb'*Mb*g;

% solve
u = A\b;

plotSol(u,50);


keyboard