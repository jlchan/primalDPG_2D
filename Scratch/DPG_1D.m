function DPG_1D

N = 2;
NT = N+2; % test order
K = 8;
r = JacobiGL(0,0,N);

% local ops
V = Vandermonde1D(N,r); 
% D = Dmatrix1D(N,r,V); M = inv(V*V'); Ks = D'*M*D;

% restriction op
rT = JacobiGL(0,0,NT); VT = Vandermonde1D(NT,rT); 
DT = Dmatrix1D(NT,rT,VT); MT = inv(VT*VT'); KT = DT'*MT*DT;
IT = Vandermonde1D(N,rT) * inv(V);

% define mesh
Np = N+1; 
VX = linspace(-1,1,K+1);
x = ones(Np,1)*VX(1:K) + 0.5*(r+1)*(VX(2:K+1)-VX(1:K));
h = diff(VX);h=h(1);

% define operators
Ks = (1/h)*kron(speye(K),KT);
M = h*kron(speye(K),MT);
RV = h*M + (1/h)*Ks;
B = (1/h)*kron(speye(K),KT*IT); % poisson

% topological ops
NpT = N*K+1; R = sparse(NpT,Np*K);
for k = 1:K
    inds1 = (k-1)*N + (1:Np);
    inds2 = (k-1)*Np + (1:Np);
    R(inds1,inds2) = speye(Np);           
end
Em = sparse(K+1,(NT+1)*K);Ep = sparse(K+1,(NT+1)*K);
for v = 2:K
    offm = (v-2)*(NT+1);
    offp = (v-1)*(NT+1);
    Em(v,offm + NT+1) = 1;      
    Ep(v,offp + 1) = 1;
end
Em(1,1) = 1;
Ep(K+1,(NT+1)*K) = 1;
Bhat = (Em-Ep)';

% stiffness
T = RV\[B*R' Bhat];
A = T'*[B*R' Bhat];
% A = R*(M+Ks)*R';

f = ones((NT+1)*K,1);     
b = [M*f; zeros(K+1,1)];
b = T'*M*f;

%bcs
A(1,:) = 0;A(:,1) = 0; A(1,1) = 1;
A(NpT,:) = 0;A(:,NpT) = 0; A(NpT,NpT) = 1;
b(1) = 0;b(NpT) = 0;
u = A\b;

plot(x(:),R'*u(1:NpT),'.')


return

