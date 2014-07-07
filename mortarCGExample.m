function err=mortarCGExample(N,Nf,mesh)

Globals2D;

% Polynomial order used for approximation 
if nargin<1
    useCG = 1;
    N = 5; % when N = even, Nf = N-1, fails?
    Nf = N-1;
    %     Read in Mesh
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
    Nv = 3;
    VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
    EToV = [3 1 2];
    K = 1;
%     [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');

else
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
end

% keyboard
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');


% Initialize solver and construct grid and metric
StartUp2D;

[M, Dx, Dy] = getBlockOps();
A = getVolOp(M,Dx,Dy);

uex = @(x,y) sin(pi*x).*sin(pi*y);
uex = project_fxn(uex,25);
f = uex*(1+2*pi^2);
% f = ones(size(x(:)));

if useCG
    [R vmapBT] = getCGRestriction();        
    A = R*A*R';
%     f = ones(Np*K,1);    
    b = R*M*f;
    
    n = size(A,1); 
    u0 = zeros(n,1); 
    % Dirichlet BC on the left
%     u0(vmapBT) = (x(vmapB) < -1+1e-7).*sqrt(1-y(vmapB).^2); 
%     vmapBT(x(vmapB) < -1+NODETOL) = []; % remove left BCs for Neumann
%     vmapBT(x(vmapB) > 1-NODETOL) = []; % remove left BCs for Neumann
    b = b - A*u0;    
    b(vmapBT) = u0(vmapBT);
    A(vmapBT,:) = 0; A(:,vmapBT) = 0;
    A(vmapBT,vmapBT) = speye(length(vmapBT));
    u = A\b;
    u = R'*u;

else
    [B, vmapBF, xf, yf, nxf, nyf] = getMortarConstraint(Nf);
                
    nU = size(B,2); % num CG nodes
    nM = size(B,1); % num mortar nodes
    O = sparse(nM,nM);     
    Am = [A B';B O];
    
    b = M*f;
    bm = [b;zeros(nM,1)];
        
    um = Am\bm;
    u = um(1:Np*K);
    
%     color_line3(x,y,u,u,'.')    
%     hold on;PlotMesh2D
end

err = u-uex;
err = sqrt(err'*M*err);

if nargin<1
    Nplot = 25;
    [xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
    Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
    xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
    yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
    figure
    color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.')
    hold on
%     plot3(x(vmapB),y(vmapB),u(vmapB),'o','markersize',8)
    plot3(x(:),y(:),u(:),'o','markersize',8)
    title(sprintf('N = %d, Nf = %d, err = %d',N, Nf, err))
end

function Vol = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;
b1 = 1; b2 = 1;ep = 1e-6;

S = -(b1*Dx+b2*Dy)'*M;
Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Vol = M + Kb + ep*Ks; % Convection-diffusion
% Vol = M+ 1e-4*Ks + Kb;  % Poisson
Vol = M + Ks;
% Vol = M + Kb;

% b1 = 1; b2 = 1;ep = 1e-1;
% L = (Dx*b1 + Dy*b2) - ep*(Dx^2 + Dy^2); 
% GLS = L'*M*L;
% tau = 1/(norm([b1 b2])*(N+1))*spdiag(sqrt(J(:)));
% As = tau*GLS;
% bs = tau*L'*M*f;

function [M, Dx, Dy] = getBlockOps()

Globals2D

blkDr = kron(speye(K),Dr);
blkDs = kron(speye(K),Ds);
blkM = kron(speye(K),MassMatrix);

M = spdiag(J(:))*blkM; % J = h^2
Dx = spdiag(rx(:))*blkDr + spdiag(sx(:))*blkDs;
Dy = spdiag(ry(:))*blkDr + spdiag(sy(:))*blkDs;


function plotNodes()
Globals2D
figure
plot(x,y,'.');hold on;
j = 1;
for k = 1:K    
    for i = 1:size(x,1)        
        off = .2*((k-1)/K);
        text(x(i,k)+off,y(i,k)+off,num2str(j),'fontsize',16)
        j = j+1;
    end
end
% PlotMesh2D

function D = spdiag(d)
n = length(d(:));
D = spdiags(d(:),0,n,n);

