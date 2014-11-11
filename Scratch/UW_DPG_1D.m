% gets 1D DPG ultra-weak operator for poisson on a uniform grid of K
% elements, with test order NT and trial order N

% UW_DPG_1D(N,NT,K)
clear

N = 4;
NT = N+2; % test order
K = 2;

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
nU = Np*K;
nF = K+1;

% define test operators
IT = kron(speye(K),IT);
DT = (1/J)*kron(speye(K),DT);
MT = J*kron(speye(K),MT);

% jumps
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

% poisson = restricted
% B = -DT'*MT*IT;
RV = J*MT + (1/J)*DT'*MT*DT;

B = [MT*IT -DT'*MT*IT;
    -DT'*MT*IT 0*IT];
RV = blkdiag(RV,RV);
Bhat = blkdiag(Bhat,Bhat);

% stiffness
T = RV\[B Bhat];
A = T'*[B Bhat];

f = ones((NT+1)*K,1);
b = T'*[zeros((NT+1)*K,1);MT*f];

%bcs - left trace
bc = 2*nU+1;
A(bc,:) = 0; A(:,bc) = 0;
A(bc,bc) = 1;
b(bc) = 0;

% right trace
bc = 2*nU + 1 + K;
A(bc,:) = 0; A(:,bc) = 0;
A(bc,bc) = 1;
b(bc) = 0;

if 0
    u = A\b;
    figure
    plot(x(:),u(Np*K+1:2*Np*K),'.-')
    xp = linspace(-1,1,100);
    hold on;plot(xp,.5*(1+xp).*(1-xp),'r-')
end

P = sparse(size(A,1),size(A,2));
for e = 1:K
    ids{e} = [(1:Np) + Np*(e-1), (1:Np) + Np*(e-1) + Np*K, 2*Np*K + (e:e+1), 2*Np*K + K + 1 + (e:e+1)]; 
    P(ids{e},ids{e}) = P(ids{e},ids{e}) + inv(A(ids{e},ids{e}));
end
cond(A)
cond(P*A)









