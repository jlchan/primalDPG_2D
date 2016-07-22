% gets 1D DPG operator for poisson on a uniform grid of K elements, wVTh
% test order NT and trial order N

function [A b ids R] = DPG_adv_1D(N,NT,K)

if nargin<3 
    N = 3;
    NT = N+2; % test order    
    K = 8;
end

r = JacobiGL(0,0,N);

% local ops
V = Vandermonde1D(N,r); 
D = Dmatrix1D(N,r,V); 
M = inv(V*V'); 

% restriction from test to trial space
rT = JacobiGL(0,0,NT); VT = Vandermonde1D(NT,rT); 
DT = Dmatrix1D(NT,rT,VT); 
MT = inv(VT*VT'); 
%KT = DT'*MT*DT;
KT = -DT'*MT;
VT = Vandermonde1D(N,rT) * inv(V); % eval basis at high order quadrature

% define mesh
Np = N+1; 
VX = linspace(-1,1,K+1);
x = ones(Np,1)*VX(1:K) + 0.5*(r+1)*(VX(2:K+1)-VX(1:K));
xT = VT*x;
h = diff(VX);h=h(1);

% define operators
Ks = (1/h)*kron(speye(K),KT);
Mtr = h*kron(speye(K),M);
M = h*kron(speye(K),MT);
RV = h*M + (1/h)*Ks;
B = (1/h)*kron(speye(K),KT*VT); % poisson

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
A(1,:) = 0; A(:,1) = 0; A(1,1) = 1; b(1) = 0;
% A(NpT,:) = 0;A(:,NpT) = 0; A(NpT,NpT) = 1; b(NpT) = 0;

if nargout==0
    u = A\b;
    figure
    plot(x(:),R'*u(1:NpT),'.-')
end

fieldIds = 1:NpT;
fluxIds = NpT+1:size(A,1);

keyboard
ids = cell(K,1);
off_field = 0;
off_flux = NpT;
for e = 1:K
    fieldIds = off_field + (1:N+1);
    fluxIds = off_flux + (1:2);
    ids{e} = [fieldIds(:); fluxIds(:)];
    off_field = off_field + N;
    off_flux = off_flux + 1;
end

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

