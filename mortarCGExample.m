function mortarCGExample

Globals2D;

% Polynomial order used for approximation 
N = 3;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');


% Initialize solver and construct grid and metric
StartUp2D;

[M, Dx, Dy] = getBlockOps();
A = getVolOp(M,Dx,Dy);

useCG = 0;
if useCG
    [R vmapBT] = getCGRestriction();        
    A = R*A*R';    
    f = ones(Np*K,1);    
    b = R*M*f;
    
    n = size(A,1); 
    u0 = zeros(n,1); 
    % Dirichlet BC on the left
%     u0(vmapBT) = (x(vmapB) < -1+1e-7).*sqrt(1-y(vmapB).^2); 
%     vmapBT(x(vmapB) < -1+NODETOL) = []; % remove left BCs for Neumann

    b = b - A*u0;    
    b(vmapBT) = u0(vmapBT);
    A(vmapBT,:) = 0; A(:,vmapBT) = 0;
    A(vmapBT,vmapBT) = speye(length(vmapBT));
    
    u = R'*(A\b);

else   
    
    [B, vmapBF, xf, yf, nxf, nyf] = getMortarConstraint(N-1);
            
    nU = size(B,2); % num CG nodes
    nM = size(B,1); % num mortar nodes
    O = sparse(nM,nM);     
    Am = [A B';B O];

    %     keyboard
    % forcing            
    %     f = zeros(size(x));
    %     f(:,round(K/2)) = 1;
    %     f = f(:);
    %     f = sin(pi*x(:)).*sin(pi*y(:));
    
    f = ones(Np*K,1);
    b = M*f;
    bm = [b;zeros(nM,1)];
    
    bmaskF = (xf < -1+NODETOL); % x = 0
    vmapBF(~bmaskF) = []; % don't impose on non-Neumann boundaries
    vmapBF = vmapBF+nU; 
    Am(vmapBF,:) = 0;Am(:,vmapBF) = 0;
    Am(vmapBF,vmapBF) = speye(length(vmapBF));
    
    um = Am\bm;
    u = um(1:Np*K);
    
%     color_line3(x,y,u,u,'.')    
%     hold on;PlotMesh2D
end

Nplot = 25;
[xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.')

function Vol = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;
b1 = 1; b2 = 1;ep = 1e-6;

S = -(b1*Dx+b2*Dy)'*M;
Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Vol = M + Kb + ep*Ks; % Convection-diffusion
Vol = Ks;  % Poisson

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

