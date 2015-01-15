function fem1D

N = 4;
K = 128;
r = JacobiGL(0,0,N);

% local ops
V = Vandermonde1D(N,r); D = Dmatrix1D(N,r,V);
M = inv(V*V'); Ks = D'*M*D;

Np = N+1;
VX = linspace(-1,1,K+1);
x = ones(Np,1)*VX(1:K) + 0.5*(r+1)*(VX(2:K+1)-VX(1:K));
h = diff(VX); h=h(1); % assume equispaced

Ks = (2/h)*kron(speye(K),Ks);
M = (h/2)*kron(speye(K),M);

NpT = N*K+1;
R = sparse(NpT,Np*K);
for k = 1:K
    inds1 = (k-1)*N + (1:Np);
    inds2 = (k-1)*Np + (1:Np);
    R(inds1,inds2) = speye(Np);
end
A = R*(Ks)*R';

f = (1+(2*pi)^2)*sin(2*pi*x(:));
b = R*M*f;

%bcs
A(1,:) = 0;A(:,1) = 0; A(1,1) = 1;
A(NpT,:) = 0;A(:,NpT) = 0; A(NpT,NpT) = 1;
b(1) = 0;b(NpT) = 0;
plot(x(:),R'*(A\b),'.')
hold on
xp = linspace(-1,1,256);
plot(xp,sin(2*pi*xp),'r-')

filename = sprintf('fem1D_N%i.mat',N);
save(filename,'A')
return
