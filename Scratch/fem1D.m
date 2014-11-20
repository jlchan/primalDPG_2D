function fem1D

N = 3;
K = 8;
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
A = R*(M + Ks)*R';

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

return

% S = inv(full(A));
% NpT = size(A,1);
% k = round(NpT/2);I1 = 1:k; I2 = k+1:NpT;
% semilogy(svd(S(I1,I1)),'.');hold on;semilogy(svd(S(I1,I2)),'ro')
% [U, flag, relres, iter, resvec] = pcg(A,b,1e-6,100,@(x) SRJ(x,A,b,NpT));
% [U, flag, relres, iter, resvec] = pcg(A,b,1e-6,200);
x = zeros(NpT,1);
[x resvec] = SRJ(x,A,b,NpT);
semilogy(resvec)

function [x rn] = SRJ(x,A,b,NpT)

% x = zeros(NpT,1);
xex = A\b;
R = A-diag(diag(A));
iD = (1./diag(A));
r = (b-A*x);
av = [50, .6];
P = [1 50];
nIter = 2000;
k = 1;
while k<nIter
    for i = 1:P(1)
        a = av(1);
        x = (1-a)*x + a*iD.*(b-R*x);
        rn(k) = norm(b-A*x);
        k = k+1;
        
        for j = 1:P(2)
            a = av(2);
            x = (1-a)*x + a*iD.*(b-R*x);
            rn(k) = norm(b-A*x);
            k = k+1;
            plot(x-xex)
            drawnow

        end
    end
end
% semilogy(rn);hold on