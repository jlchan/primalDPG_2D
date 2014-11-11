% gets 1D DPG ultra-weak operator for poisson on a uniform grid of K
% elements, with test order NT and trial order N

% UW_DPG_1D(N,NT,K)

N = 2;
NT = N+2; % test order
K = 8;

r = JacobiGL(0,0,N);

% local ops
V = Vandermonde1D(N,r);
D = Dmatrix1D(N,r,V); M = inv(V*V'); Ks = D'*M*D;

% restriction from test to trial space
rT = JacobiGL(0,0,NT); VT = Vandermonde1D(NT,rT);
DT = Dmatrix1D(NT,rT,VT); MT = inv(VT*VT'); ST = -DT'*MT;
IT = Vandermonde1D(N,rT) * inv(V);

% define mesh
Np = N+1;
VX = linspace(-1,1,K+1);
x = ones(Np,1)*VX(1:K) + 0.5*(r+1)*(VX(2:K+1)-VX(1:K));
xT = ones(NT+1,1)*VX(1:K) + 0.5*(rT+1)*(VX(2:K+1)-VX(1:K));
h = diff(VX);h=h(1); J = h/2;

% define test operators
ST = kron(speye(K),-DT'*MT);
DT = (1/J)*kron(speye(K),DT);
MT = J*kron(speye(K),MT);
RV = J*MT + (1/J)*DT'*MT*DT;

% poisson = restricted
B = ST*kron(speye(K),IT); 

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
T = RV\[B Bhat];
A = T'*[B Bhat];

f = ones((NT+1)*K,1);
b = [MT*f; zeros(K+1,1)];
b = T'*MT*f;

%bcs
nU = Np*K;
nF = K+1;
bc = nU+1;
A(bc,:) = 0; A(:,bc) = 0;
A(bc,bc) = 1;
% A(NpT,:) = 0;A(:,NpT) = 0; A(NpT,NpT) = 1;
b(1) = 0;
%b(NpT) = 0;

% fieldInds = 1:NpT;
% fluxInds = NpT+1:size(A,1);

u = A\b;
figure
plot(x(:),u(1:Np*K),'.-')


