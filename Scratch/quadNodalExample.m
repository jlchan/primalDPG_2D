function quadNodalExample

% set up 1D operations
N = 10; Np = N+1; Np2 = Np^2;
[r] = JacobiGL(0,0,N);
[r s] = meshgrid(r);r = r(:);s = s(:);

% get derivative matrices
V = Vander2D(r,s,N); invV = inv(V); [Dx Dy] = Grad2D(r,s,invV,N);

% define mass/stiffness matrices
M = invV'*invV;
K = Dx'*M*Dx + Dy'*M*Dy;
S = (Dx+.5*Dy)'*M*(Dx+.5*Dy);

f = r.*(1-s);
f = zeros(Np2,1);
f(round(Np2/2)) = 1;
A = M + .0*K + S;
b = M*f; % f = 1
% bcs
tol = 1e-14;
onB = @(x,y) 0*(abs(1-x.^2) < tol) | (abs(1-y.^2) < tol);
bMap = find(onB(r,s));
for i = 1:length(bMap)
%     A(bMap(i),:) = A(bMap(i),:)*0; A(:,bMap(i)) = A(:,bMap(i))*0;
%     A(bMap(i),bMap(i)) = 1; b(bMap(i)) = 0;
end

% solve
u = A\b;

% get smallest generalized eigval for coercivity
gam = min(eig(A,M));

% uniform plotting grid
[ru] = linspace(-1,1,120); [ru su] = meshgrid(ru);ru = ru(:);su = su(:);
Vu = Vander2D(ru,su,N); Vp = Vu*invV; 
color_line3(ru,su,Vp*u,Vp*u,'.');
hold on
color_line3(r(bMap),s(bMap),u(bMap),u(bMap),'o');
title(['Coercivity const = ', num2str(gam)])
figure;
color_line3(ru,su,Vp*f,Vp*f,'.')

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
