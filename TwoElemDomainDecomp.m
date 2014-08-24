% function to test domain decomposition methods using global spectral
% methods on each subdomain 

function TwoElemDomainDecomp

% set up 1D operations
N = 12; Np = N+1; Np2 = Np^2;
Nf = 2; Nfp = Nf+1;

[r1D] = JacobiGL(0,0,N);
[r s] = meshgrid(r1D);r = r(:);s = s(:);

% get derivative/mass matrices
V = Vander2D(r,s,N); invV = inv(V); 
[Dx Dy] = Grad2D(r,s,invV,N);
M = invV'*invV;

% define subdomain stiffness matrix/forcing
% A = M + Dx'*M*Dx + Dy'*M*Dy;
% L = (Dx + .5*Dy);
% A = M + L'*M*L;
A = M + Dx'*M*Dx + Dy'*M*Dy;

tol = 1e-14;
V1D = Vandermonde1D(N,r1D); invV1D = inv(V1D); M1D = invV1D'*invV1D;

% interface conditions
u1I = find(abs(1-r) < tol);
u2I = find(abs(1+r) < tol) + Np2; % corresponding to 2nd block 
E1 = sparse(1:length(u1I),u1I,1,Np,2*Np2);
E2 = sparse(1:length(u2I),u2I,1,Np,2*Np2);
Ef = E1 - E2; %normal jump in u over interface

% interface condition
if (Nf==0)
    Rf = ones(Np,1); % interpolate constant to multiple nodes
else
    r1Df = JacobiGL(0,0,Nf); 
    V1Df = Vandermonde1D(Nf,r1Df); V1Dfx = Vandermonde1D(Nf,r1D);  
    Rf = V1Dfx*inv(V1Df);    
end
B = Rf'*M1D*Ef;

% keyboard
blkA = blkdiag(A,A);
g = Np2;

blkA = blkA + g*(E1'*(M1D*E1) - E2'*(M1D*E2));
blkM = blkdiag(M,M);
% blkDx = blkdiag(Dx,Dx);
% blkDy = blkdiag(Dy,Dy);
A = [blkA B'
    B zeros(Nfp)];

x = [r; r+2];
y = [s; s];
a = 1;
uex = @(x,y) sin(a*pi*x).*cos(pi*y);
f = (1+(1+a^2)*pi^2)*uex(x,y);
% f = [-(s<0); (s>0)];
% f = [-(r<0); (r>0)];
% f = ones(2*Np2,1);
b = blkM*f; 
bf = zeros(Nfp,1); %bf(round(Nfp/2)) = 1;
b = [b; bf];

% homogeneous bcs
u1BCs = find(r < -1+tol);
u2BCs = find(r > 1-tol);
bMap = [u1BCs; u2BCs + Np2];
A(bMap,:) = 0; A(:,bMap) = 0;
A(bMap,bMap) = speye(length(bMap)); 
b(bMap) = 0;

% solve
u = A\b;
u1 = u(1:Np2);
u2 = u(Np2+1:2*Np2);
u = u(1:Np2*2);

% uniform plotting grid
[ru1D] = linspace(-1,1,120); [ru su] = meshgrid(ru1D);ru = ru(:);su = su(:);
Vu = Vander2D(ru,su,N); Iu = Vu*invV; 
figure
subplot(2,1,1)
color_line3(ru,su,Iu*u1,Iu*u1,'.');
color_line3(ru+2+.075,su,Iu*u2,Iu*u2,'o');

[rcub scub wcub] = cub2D(25);
Vc = Vander2D(rcub,scub,N); Ic = Vc*invV;
e1 = uex(rcub,scub) - Ic*u1;
err1 = sqrt(e1'*(wcub.*e1));
e2 = uex(rcub+2,scub) - Ic*u2;
err2 = sqrt(e2'*(wcub.*e2));

e1 = uex(r,s) - u1;
e2 = uex(r+2,s) - u2;
% color_line3(ru,su,Iu*e1,Iu*e1,'.');
% color_line3(ru+2+.075,su,Iu*e2,Iu*e2,'o');


B = A(2*Np2+(1:Nfp),1:2*Np2);
A = A(1:2*Np2,1:2*Np2);
S = full(B*(A\B'));
cA = condest(A)
cS = cond(S)
title(['Norm of err 1 = ' num2str(err1), ', norm of err 2 = ' num2str(err2), ', logcond of blk A = ', num2str(log(cA)/log(10))])


Vu = Vandermonde1D(N,ru1D);Iu1D = Vu*invV1D;
subplot(2,1,2)
ju = Ef*u;
plot(ru1D,Iu1D*ju);
title(['Norm of jump at interface = ' num2str(sqrt(ju'*M1D*ju))])


function V = Vander2D(x,y,N)

Np = N+1;
V = zeros(length(x),Np^2);
k = 1;

for i = 0:N % over x
    Px = JacobiP(x,0,0,i);
    for j = 0:N % over y        
        V(:,k) = Px.*JacobiP(y,0,0,j);
        k = k+1;
    end
end

function [Dx Dy] = Grad2D(x,y,invV,N)

Np = N+1;
Dx = zeros(length(x),Np^2);
Dy = zeros(length(x),Np^2);
k = 1;
for i = 0:N % over x
    Px = JacobiP(x,0,0,i);
    DPx = GradJacobiP(x,0,0,i);
    for j = 0:N % over y        
        Dx(:,k) = DPx.*JacobiP(y,0,0,j);
        Dy(:,k) = Px.*GradJacobiP(y,0,0,j);
        k = k+1;
    end
end
Dx = Dx*invV;
Dy = Dy*invV;

% interp cubature 
function [rcub scub cubW] = cub2D(Ncub)

[x,w] = JacobiGQ(0,0,Ncub);
[rcub scub] = meshgrid(x);
[wr ws] = meshgrid(w);

rcub = rcub(:);
scub = scub(:);
cubW = wr(:).*ws(:);



