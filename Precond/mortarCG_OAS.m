function err=mortarCGExample(N,Nf,mesh)

Globals2D;

% Polynomial order used for approximation 
if nargin<1
    useCG = 0;
    N = 6; % when N = even, Nf = N-1, fails?
    Nf = 3;
    %     Read in Mesh
    mesh = 'squarereg.neu';
    mesh = 'Maxwell025.neu';
%     mesh = 'Maxwell1.neu';
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
    %     Nv = 3; VX = VX(EToV(1,:)); VY = VY(EToV(1,:)); EToV = [3 1 2]; K = 1;
    
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
    
    % build OAS preconditioner (1 elem overlap)
    [r c] = find(R); r = reshape(r,Np,K); c = reshape(c,Np,K);
    Ak = cell(K,1);Aki = cell(K,1);
    for k = 1:K
        nbr = k;%unique([k EToE(k,:)]);
        inds = unique(r(:,nbr));
        Aki{k} = inds;
        Ak{k} = eye(numel(inds));%A(inds,inds);                   
%         clf
%         plotVerts;hold on
%         ci = unique(c(:,nbr));
%         plot(x(ci),y(ci),'ro')
%         title(['k = ', num2str(k)])
%         pause
    end    
%     u = OAS(b,Ak,Aki);
%     plotSol(R'*u,25)    
%     keyboard

    % coarse grid solver
    Rs = diag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
    [Rc1 Ir1 vmapBT1 xr1 yr1] = pRestrictCG(N,1); % interp down
    Rp1 = Rc1*Ir1';
    R1 = Rs*Rp1'; % interp down to P1
    [r c] = find(R1);  bci1 = ismember(r,vmapBT);
    vmapBT1 = unique(c(bci1));
    % build boundary conditions into R1
    %     R1(vmapBT,:) = 0;  R1(:,vmapBT1) = 0;
    %     R1(vmapBT1,vmapBT1) = speye(length(vmapBT1));
    A1 = R1'*A*R1; A1(vmapBT1,:) = 0;A1(:,vmapBT1) = 0;A1(vmapBT1,vmapBT1) = speye(numel(vmapBT1));
    b1 = R1'*b; b1(vmapBT1) = 0; % FIGURE OUT HOW TO USE u0!!!
    P1 = @(b) R1*((R1'*A*R1)\(R1'*b));    
    
    % build coarse grid solver
    [u, flag, relres, iter, resvec] = pcg(A,b,1e-6,50,@(x) OAS(x,Ak,Aki) + P1(x));
    
    semilogy(resvec,'.-'); hold on;
    title('OAS for CG')
    keyboard
    u = R'*u;    
else
%     [B, vmapBF, xfb, yfb, nxf, nyf fmap xf yf] = getMortarConstraint(Nf);    
    [B vmapBF xf yf nxf nyf fmap fpairs] = getMortarConstraint(Nf);        
    
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
        
    % build OAS preconditioner
    S = B*(A\B'); bS = B*(A\b);    
    FToE = getFToE(fpairs);
    NfacesU = size(fpairs,2);
    
    fmapU = reshape(1:(Nf+1)*NfacesU,Nf+1,NfacesU);            
%         M = sparse(size(S,1),size(S,2));
    Sf = cell(NfacesU,1);
    Sfi = cell(NfacesU,1);
    for f = 1:NfacesU
%         inds = unique(fmap(:,FToE(f,:)));
        inds = fmapU(:,f); % no overlap
%                 M(inds,inds) = M(inds,inds) + inv(S(inds,inds));
        Sf{f} = S(inds,inds);
        Sfi{f} = inds;
    end    
    
    % build coarse grid solver
    If1 = ones(Nf+1,1)/sqrt(Nf+1); % interp from constant to Nf+1
    If1 = kron(speye(NfacesU),If1);    
    P1 = @(bS) If1*((If1'*S*If1)\(If1'*bS));
    
    [f, flag, relres, iter, resvec] = pcg(S,bS,1e-6,50,@(x) OAS(x,Sf,Sfi) + P1(x));
    semilogy(resvec,'.-');    hold on
    title(sprintf('OAS for mortars with N = %d, mesh = %s',N,mesh))    
    
    u = A\(b-B'*f);
end

err = u-uex;
err = sqrt(err'*M*err);

if nargin<1
    plotSol(u,25);    
    title(sprintf('N = %d, Nf = %d, err = %d',N, Nf, err))        
end
keyboard

function Vol = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;
% b1 = 1; b2 = 1;ep = 1e-6;

% S = -(b1*Dx+b2*Dy)'*M;
% Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Vol = M + Kb + ep*Ks; % Convection-diffusion
% Vol = M+ 1e-4*Ks + Kb;  % Poisson
Vol = M + Ks;
% Vol = M + Kb;


function b = OAS(x,Af,Afi)
d = zeros(size(x));
for i = 1:length(Afi)
    d(Afi{i}) = d(Afi{i}) + 1;
end
d = 1./d; % scale
d = d.^0;

b = zeros(size(x));
for i = 1:length(Afi)    
    b(Afi{i}) = b(Afi{i}) +  d(Afi{i}).*(Af{i}\x(Afi{i}));    
end

