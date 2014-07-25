% experiments w/tensor products

Globals2D; FaceGlobals2D

N = 2; % when N = even, Nf = N-1, fails?
Nf = N-2; % = N trial
Nt = N-1; % = 

[Nv, VX, VY, K, EToV] = QuadMesh2D(1);

% Initialize solver and construct grid and metric
StartUp2D; FaceStartUp2D

Dr = full(Dr);Ds = full(Ds);
[M, Dx, Dy] = getBlockOps();
Np1D = N+1;
r1D = x(1:Np1D);
V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D');
D1D = Dmatrix1D(N,r1D,V1D);

DxK = kron(eye(Np1D),D1D);
DyK = kron(D1D,eye(Np1D));
norm(full(Dx)-DxK,'fro')
norm(full(Dy)-DyK,'fro')

K1D = D1D'*M1D*D1D;
KxK = kron(M1D,K1D) % x is stored in the second half
Kx = Dx'*MassMatrix*Dx


