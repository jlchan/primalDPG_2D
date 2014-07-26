% tensor product fast inversion 

Globals2D
[Nv, VX, VY, K, EToV] = QuadMesh2D(1);
N = 4;
StartUp2D

Dr = full(Dr);Ds = full(Ds);
M = MassMatrix;
f = ones(Np,1);
b = M*f;

Np1D = N+1;
r1D = x(1:Np1D);
V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D');
D1D = Dmatrix1D(N,r1D,V1D);
invM = V1D*V1D';

K1D = D1D'*M1D*D1D;

% 1D BCs
K1D(1,:) = 0;K1D(:,1) = 0;K1D(1,1) = 1; 
K1D(Np1D,:) = 0;K1D(:,Np1D) = 0;K1D(Np1D,Np1D) = 1;
M1D(1,:) = 0;M1D(:,1) = 0;M1D(1,1) = 1;
M1D(Np1D,:) = 0;M1D(:,Np1D) = 0;M1D(Np1D,Np1D) = 1;

% b(vmapB) = 0;
b = reshape(b,Np1D,Np1D);
bM = invM*b*invM;
bM(vmapB) = 0;
I = eye(Np1D);
Kh = kron(I,K1D/M1D) + kron(M1D\K1D,I);
plotSol(Kh\bM(:),50)

[Q L] = eig(K1D,M1D); % they're eigenvectors for both inv(M)*K and K*inv(M)
iQ = Q';
% [P ~] = eig(K1D*inv(M1D)); 
P = M1D*Q;
iP = P';

[Lam LamT] = meshgrid(diag(L));
iLam = 1./(Lam+LamT);

u = P*(iLam.*(iP*bM*iQ))*Q;
u = reshape(u,Np1D^2,1);
plotSol(u,50);
