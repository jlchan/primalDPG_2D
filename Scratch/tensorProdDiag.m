% tensor product fast inversion 

Globals2D
[Nv, VX, VY, K, EToV] = QuadMesh2D(1);
N = 10;
StartUp2D

Dr = full(Dr);Ds = full(Ds);
M = MassMatrix;
f = ones(Np,1);
% f = randn(Np,1);
b = M*f;
b(vmapB)=0;

%% ref sol
Kh = Dr'*M*Dr + Ds'*M*Ds;
Kh(vmapB,:) = 0;Kh(:,vmapB) = 0;Kh(vmapB,vmapB) = speye(length(vmapB));
uR = Kh\b;

%% exact tensor product sol
Np1D = N+1;
r1D = JacobiGL(0,0,N);
V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D');
D1D = Dmatrix1D(N,r1D,V1D);
K1D = D1D'*M1D*D1D;

% 1D BCs8
K1D(1,:) = 0;K1D(:,1) = 0;K1D(1,1) = 1; 
M1D(1,:) = 0;M1D(:,1) = 0;M1D(1,1) = 1;
K1D(Np1D,:) = 0;K1D(:,Np1D) = 0;K1D(Np1D,Np1D) = 1;
M1D(Np1D,:) = 0;M1D(:,Np1D) = 0;M1D(Np1D,Np1D) = 1;

a = 1;
Kk = a*kron(M1D,M1D) + kron(M1D,K1D) + kron(K1D,M1D);
uT=Kk\b;

%% diagonalization sol
% only need one set of eigvectors for isotropic poisson
% kron(M,K) + kron(K,M) reversible

[Q, L] = eig(K1D,M1D); % they're eigenvectors for both inv(M)*K and K*inv(M)
Q = Q*diag(1./sqrt(diag(Q'*M1D*Q))); %normalize in M
[Lam LamT] = meshgrid(diag(L));
iLam = 1./(a + Lam + LamT);

% % equiv to these 1D dense ops
% kronL = kron(.5*a*eye(Np1D) + L,eye(Np1D)) + kron(eye(Np1D),.5*a*eye(Np1D)+L);
% u = kron(Q,Q)*( kronL\(kron(Q',Q')*b(:)));

b = reshape(b,Np1D,Np1D);
u = Q*(iLam.*(Q'*b*Q))*Q';
uD = u(:);

plotSol(uD,75);
plotSol(uT-uD,75)
title(['norm of diff = ', num2str(norm(uT-uD))])