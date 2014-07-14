function err=mortarCGExample(N,Nf,mesh)

Globals2D;

% Polynomial order used for approximation 
if nargin<1
    useCG = 0;
    N = 3; % when N = even, Nf = N-1, fails?
    Nf = 1;
    %     Read in Mesh
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
%     Nv = 3;
%     VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
%     EToV = [3 1 2];
%     K = 1;
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
%     [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');

else
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
end

% Initialize solver and construct grid and metric
StartUp2D;

[M, Dx, Dy] = getBlockOps();
A = getVolOp(M,Dx,Dy);

uex = @(x,y) sin(pi*x).*sin(pi*y);
uex = project_fxn(uex,25);
f = uex*(1+2*pi^2);
f = ones(size(x(:)));

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
%     [B, vmapBF, xfb, yfb, nxf, nyf fmap xf yf] = getMortarConstraint(Nf);    
    [B vmapBF xf yf nxf nyf fpairs] = getMortarConstraint(Nf);        
    
    xfb = xf(vmapBF);yfb = yf(vmapBF);
    nxf = nxf(vmapBF);nyf = nyf(vmapBF);        
    nU = size(B,2); % num CG nodes
    nM = size(B,1); % num mortar nodes
    O = sparse(nM,nM);     
    Am = [A B';B O];
    
    b = M*f;
    bm = [b;zeros(nM,1)];
        
    um = Am\bm;
    u = um(1:Np*K);
    f = um(Np*K+1:end);        
end

err = u-uex;
err = sqrt(err'*M*err);

if nargin<1
    plotSol(u,25);    
    plotSol(sqrt((Dx*u).^2 + (Dy*u).^2),25)    
    title(sprintf('N = %d, Nf = %d, err = %d',N, Nf, err))        
end
keyboard

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
