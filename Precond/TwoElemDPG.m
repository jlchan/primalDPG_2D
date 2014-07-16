% function to test domain decomposition methods using global spectral
% methods on each subdomain 

function TwoElemDPG
tol = 1e-14;

NTr = 10; NpTr = NTr+1; NpTr2 = NpTr^2;
N = NTr+1; Np = N+1; Np2 = Np^2; % test order/dofs 
Nf = NTr; Nfp = Nf+1;

% test 
[r1D] = JacobiGL(0,0,N);
[r s] = meshgrid(r1D);r = r(:);s = s(:);
V = Vander2D(r,s,N); invV = inv(V); 
M = invV'*invV; [Dx Dy] = Grad2D(r,s,invV,N); 
ATest = M + Dx'*M*Dx + Dy'*M*Dy;

% trial restriction
[r1DTr] = JacobiGL(0,0,NTr);
[rTr sTr] = meshgrid(r1DTr);rTr = rTr(:);sTr = sTr(:);
VTr = Vander2D(rTr,sTr,NTr); invVTr = inv(VTr);
VTr2Test = Vander2D(r,s,NTr);
ITr = VTr2Test*invVTr; % trial to enriched space interpolant
% A = Dx'*M*Dx + Dy'*M*Dy;
A  = M*Dx + .1*(Dx'*M*Dx + Dy'*M*Dy);
A = A*ITr;

% interface flux - eval @ test order
V1D = Vandermonde1D(N,r1D); invV1D = inv(V1D); M1D = invV1D'*invV1D;
v1I = find(abs(1-r) < tol); v2I = find(abs(1+r) < tol) + Np2; % corresponding to 2nd block 
E1 = sparse(1:length(v1I),v1I,1,Np,2*Np2);
E2 = sparse(1:length(v2I),v2I,1,Np,2*Np2);
Ef = E1 - E2; %normal jump in v over interface
if (Nf==0)
    Rf = ones(Np,1); % interpolate constant to multiple nodes
else
    r1Df = JacobiGL(0,0,Nf); 
    V1Df = Vandermonde1D(Nf,r1Df); V1Dfx = Vandermonde1D(Nf,r1D);  
    Rf = V1Dfx*inv(V1Df);    
end
Bh = (Rf'*M1D*Ef)';

% restrict to CG space
NpTrK = NpTr2*2 - NpTr; % total C0 dofs
u1I = find(abs(1-rTr) < tol); 
u2I = find(abs(1+rTr) < tol) + NpTr2; % corresponding to 2nd block 
u1Int = find(rTr<1); u2Int = find(rTr>-1) + NpTr2; % interior nodes
vmap = zeros(2*NpTr2,1);  
vmap(u1Int) = 1:numel(u1Int);
vmap(u2Int) = numel(u1Int) + (1:numel(u2Int));
vmap(u1I) = numel(u1Int)+numel(u2Int) + (1:NpTr);
vmap(u2I) = numel(u1Int)+numel(u2Int) + (1:NpTr);
R = sparse(vmap,1:NpTr2*2,1);

% build lsq problem
x = [r; r+2];
y = [s; s];
a = 1;
uex = @(x,y) sin(a*pi*x).*cos(pi*y);
f = (1+(1+a^2)*pi^2)*uex(x,y);
f = ones(2*Np2,1);
f = x.*y+2;

blkM = blkdiag(M,M);
blkA = blkdiag(A,A)*R';
blkATest = blkdiag(ATest,ATest);
bMap = find(x < -1+tol | x > 3-tol | abs(1-y.^2)<tol); % test function BCs - isolate only middle flux
blkATest(bMap,:) = 0;blkATest(:,bMap) = 0;
blkATest(bMap,bMap) = speye(length(bMap));
A = [blkA Bh];
T = blkATest\A;
A = T'*A;
b = T'*blkM*f; 

% homogeneous bcs
xTr = [rTr; rTr + 2]; 
yTr = [sTr; sTr];
xTr = xTr([u1Int; u2Int; u1I]); yTr = yTr([u1Int; u2Int; u1I]);
bMap = find(xTr < -1+tol | xTr > 3-tol | abs(1-yTr.^2)<tol);
A(bMap,:) = 0; A(:,bMap) = 0;
A(bMap,bMap) = speye(length(bMap)); 
b(bMap) = 0;

useDD = 1;
if useDD
    figure
    Reduce = diag(1./sum(R,2))*R;
    DxTr = ITr'*Dx*ITr;
    uh = rand(Nfp,1);
    for k = 1:25
        I1 = find(xTr<1+tol);
        I2 = find(xTr>1-tol);
        A1 = A(I1,I1);
        A2 = A(I2,I2);
        
        I3 = size(A,1)-Nfp+1:size(A,1);
        B1 = A(I1,I3);
        B2 = A(I2,I3);
        C = A(I3,I3);
        b1 = b(I1)-B1*uh;
        b2 = b(I2)-B2*uh;
        u1 = A1\b1;
        u2 = A2\b2;
        A1 = A(1:NpTrK,1:NpTrK);
        B = A(1:NpTrK,NpTrK + (1:Nfp));
        b1 = b(1:NpTrK);
        u = A1\(b1 - B*uh);
        b3 = b(I3);
        uh = C\(b3 - B'*u);        
        
        % mix
%         iface = I1(ismember(I1,I2));
%         dudx = Reduce*[DxTr*u1; -0*DxTr*u2];
%         a = 1;
%         uh = a*uh + (1-a)*(-dudx(iface));

%         uh = C\(b3 - B1'*u1 - B2'*u2);       


        % test reduce + prlong        
%         u = Reduce*[u1;u2];
        u = R'*u(1:NpTrK);
        u1 = u(1:NpTr2);
        u2 = u(NpTr2+1:2*NpTr2);
        [ru1D] = linspace(-1,1,50); [ru su] = meshgrid(ru1D);ru = ru(:);su = su(:);
        Vu = Vander2D(ru,su,NTr); Iu = Vu*invVTr;
        
        clf
        subplot(2,1,1)
        color_line3(ru,su,Iu*u1,Iu*u1,'.');hold on
        color_line3(ru+2,su,Iu*u2,Iu*u2,'.');
        title(['iter k = ',num2str(k)])
        [r1Df] = JacobiGL(0,0,Nf);
        V1Df = Vandermonde1D(Nf,r1Df); invV1Df = inv(V1Df);
        Vu = Vandermonde1D(Nf,ru1D);
        Iu1D = Vu*invV1Df;
        subplot(2,1,2)
        plot(ru1D,Iu1D*uh);
        title('Interface flux')
        
        pause
    end
end

% solve
u = A\b;
uh = u(NpTrK + (1:Nfp));
u = R'*u(1:NpTrK);
u1 = u(1:NpTr2);
u2 = u(NpTr2+1:2*NpTr2);

% uniform plotting grid
[ru1D] = linspace(-1,1,50); [ru su] = meshgrid(ru1D);ru = ru(:);su = su(:);
Vu = Vander2D(ru,su,NTr); Iu = Vu*invVTr; 
figure
subplot(2,1,1)
color_line3(ru,su,Iu*u1,Iu*u1,'.');hold on
color_line3(ru+2,su,Iu*u2,Iu*u2,'.');

% plot interface flux
[r1Df] = JacobiGL(0,0,Nf);
V1Df = Vandermonde1D(Nf,r1Df); invV1Df = inv(V1Df);
Vu = Vandermonde1D(Nf,ru1D);
Iu1D = Vu*invV1Df;
subplot(2,1,2)
plot(ru1D,Iu1D*uh);
title('Interface flux')


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



