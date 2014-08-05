% gets 1D DPG operator for poisson on a uniform grid of K elements, with
% test order NT and trial order N

function [A Mtr] = DPG_1D(N,NT,K)

if nargin<3
    N = 2;
    NT = N+2; % test order    
    K = 8;
end

r = JacobiGL(0,0,N);

% local ops
V = Vandermonde1D(N,r); 
D = Dmatrix1D(N,r,V); M = inv(V*V'); Ks = D'*M*D;

% restriction from test to trial space
rT = JacobiGL(0,0,NT); VT = Vandermonde1D(NT,rT); 
DT = Dmatrix1D(NT,rT,VT); MT = inv(VT*VT'); KT = DT'*MT*DT;
IT = Vandermonde1D(N,rT) * inv(V);

% define mesh
Np = N+1; 
VX = linspace(-1,1,K+1);
x = ones(Np,1)*VX(1:K) + 0.5*(r+1)*(VX(2:K+1)-VX(1:K));
xT = ones(NT+1,1)*VX(1:K) + 0.5*(rT+1)*(VX(2:K+1)-VX(1:K));
h = diff(VX);h=h(1);

% define operators
Ks = (1/h)*kron(speye(K),KT);
Mtr = h*kron(speye(K),M);
M = h*kron(speye(K),MT);
RV = h*M + (1/h)*Ks;
B = (1/h)*kron(speye(K),KT*IT); % poisson

% topological operators
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

f = ones((NT+1)*K,1);
b = [M*f; zeros(K+1,1)];
b = T'*M*f;

%bcs
A(1,:) = 0;A(:,1) = 0; A(1,1) = 1;
A(NpT,:) = 0;A(:,NpT) = 0; A(NpT,NpT) = 1;
b(1) = 0;b(NpT) = 0;

if nargin<3
    u = A\b;
    figure
    plot(x(:),R'*u(1:NpT),'.-')
end

fieldInds = 1:NpT;
fluxInds = NpT+1:size(A,1);
% [X Y] = meshgrid(x(:));
% % plot(X,Y,'.')
% 
% K2 = kron(A,M) + kron(M,A);
% % T2 = (kron(T'*M,M) + kron(M,T'*M));
% 
% b2 = zeros(size(K2,1),1);
% b2(1)= 1;
% U = K2\b2;
% u = U(1:NpT^2);
% color_line3(X(:),Y(:),u,u,'.');

return

