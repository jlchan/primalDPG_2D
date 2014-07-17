function err= mortarCG_OAS(N,Nf,mesh)
%function mortarCG_OAS

Globals2D;
FaceGlobals2D

useCG = 0;

% Polynomial order used for approximation
if nargin<3
    N = 8;  % when N = even, Nf = N-1, fails?
    Nf = 6;
    %     Read in Mesh
    mesh = 'squarereg.neu';
    mesh = 'Maxwell05.neu';
%         mesh = 'Maxwell1.neu';
    % mesh = 'lshape.neu';
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
    %Nv = 3; VX = VX(EToV(1,:)); VY = VY(EToV(1,:)); EToV = [3 1 2]; K = 1;
end
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);


% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

[M, Dx, Dy] = getBlockOps();
AK = getVolOp(M,Dx,Dy);

uex = @(x,y) sin(pi*x).*sin(pi*y);
uex = project_fxn(uex,25);
f = uex*(1+2*pi^2);
f = ones(size(x(:)));

if useCG
    [R vmapBT] = getCGRestriction();
    A = R*AK*R';
    b = R*M*f;
    
    n = size(A,1);
    u0 = zeros(n,1);
    % Dirichlet BC on the left
    u0(vmapBT) = (x(vmapB) < -1+1e-7).*sqrt(1-y(vmapB).^2);
    %     vmapBT(x(vmapB) < -1+NODETOL) = []; % remove left BCs for Neumann
    b = b - A*u0;
    b(vmapBT) = u0(vmapBT);
    A(vmapBT,:) = 0; A(:,vmapBT) = 0;
    A(vmapBT,vmapBT) = speye(length(vmapBT));
    %     u = A\b;
    
    %     levels = agmg_setup(A1);
    %     [x flag relres iter resvec1] = agmg_solve(levels, b1, 50, 1e-6);
    %     semilogy(resvec1,'.-')
    
    Pre = buildOAS_CG(R,A);
    [u, flag, relres, iter, resvec] = pcg(A,b,1e-7,75,@(x) Pre(x));
    %[u, flag, relres, iter, resvec] = pcg(A,b,1e-6,75,@(x) OAS(x,Ak,Aki)+ P1(x));
    %     [u, flag, relres, iter, resvec] = gmres(A,b,[],1e-6,50,@(x) OAS(x,Ak,Aki) + P1(x));
    
    semilogy(resvec,'ro-'); hold on;
    title('OAS for CG')
    
    u = R'*u;
else
    B = getMortarConstraint();
    
    xfb = xf(fmapB);yfb = yf(fmapB);
    nxf = nxf(fmapB);nyf = nyf(fmapB);
    nU = Np*K; % num CG nodes
    nM = Nfrp*NfacesU; % num mortar nodes
    O = sparse(nM,nM);
    Am = [AK B';B O];
    b = M*f;
    bm = [b;zeros(nM,1)];
    
    % homogeneous BCs on V are implied by mortars.
    % BCs on mortars removes BCs on test functions.    
    bmaskf = (xfb < -1 + NODETOL) & (nxf < -NODETOL);
    fmapB(~bmaskf) = [];         

    bci = nU + fmapB; % skip over u dofs
    bm(bci) = 0;
    Am(bci,:) = 0; Am(:,bci)=0;
    Am(bci,bci) = speye(length(bci));

%     U = Am\bm;
%     u = U(1:nU);
    
    C = Am(nU + (1:nM),nU + (1:nM));
    B = Am(nU + (1:nM),1:nU);
    AK = Am(1:nU,1:nU);
    b = bm(1:nU);
    c = bm(nU+(1:nM));

    %     bm = [b;zeros(nM,1)];
    %     um = Am\bm;
    %     u = um(1:Np*K);
    %     f = um(Np*K+1:end);
    
    S = C-B*(AK\B'); bS = c-B*(AK\b);

    Pre = buildOAS_mortar(S);
    %     levels1 = agmg_setup(S1);
    %     [x flag relres iter resvec1] = agmg_solve(levels1, bS1, 50, 1e-6);
    %     semilogy(resvec1,'.-')
    %     keyboard
    
    [f, flag, relres, iter, resvec] = pcg(-S,-bS,1e-6,75,@(x) Pre(x));
    
    %[f, flag, relres, iter, resvec] = pcg(S,bS,1e-6,75,@(x) OAS(x,Sf,Sfi) + P1(x));
    %     [f, flag, relres, iter, resvec] = gmres(S,bS,[],1e-6,50,@(x) OAS(x,Sf,Sfi) + P1(x));
    semilogy(resvec,'r.-');    hold on
    title(sprintf('OAS for mortars with N = %d, mesh = %s',N,mesh))

    u = AK\(b-B'*f);
end

err = u-uex;
err = sqrt(err'*M*err);
% return
if nargin<1
    plotSol(u,25);
    title(sprintf('N = %d, Nf = %d, err = %d',N, Nf, err))
end
% keyboard

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

b = zeros(size(x));
for i = 1:length(Afi)
    b(Afi{i}) = b(Afi{i}) +  d(Afi{i}).*(Af{i}\x(Afi{i}));
end

