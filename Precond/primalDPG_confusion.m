function [A b nU nM Np Rp Irp, M, fpairs] = primalDPG_confusion(mesh,Ntrial,Ntest,Nflux,plotFlag,b,epsilon,quads)

if nargin<8
    noJacobians = 0;
end
Globals2D
FaceGlobals2D

if nargin<8
    if quads
        [Nv, VX, VY, K, EToV] = QuadMesh2D(8);
    else
        [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
        Nv = 3; VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
        EToV = [3 1 2];    K = 1;
        [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell00625.neu');        
    end
    Ntrial = 2;
    N = Ntrial+2; % Ntest
    Nf = Ntrial;
    plotFlag = 1;
    b = 1; epsilon = 1e-4; % poisson, w/beta = 0
else
    N = Ntest;
    Nf = Nflux;

    % Read in Mesh
    if quads
        [Nv, VX, VY, K, EToV] = QuadMesh2D(2^(1+mesh)); % assume mesh = 1 at smallest - 4, 8, 16
    else
        [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
    end    
end

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

global b1
global b2
global ep
b1 = b; b2 = 0; ep = epsilon;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = ones(Np*K,1);

% f = y(:)<=0;
% f = sin(pi*x(:)).*sin(pi*y(:));

[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xr yr] = pRestrictCG(N,Ntrial); % restrict test to trial space
Rr = Rp*Irp';
% Rr = Irp'; warning('discontinuous discretization!')
Bhat = getMortarConstraint();
xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);nyfb = nyf(fmapB);

B = BK*Rr';   % form rectangular bilinear form matrix
% B = BK;

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
if 1
    tic
    for i = 1:K % independently invert
        inds = (i-1)*Np + (1:Np);
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
        disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
    disp(['time for test function computation = ', num2str(toc)])
else
    disp('parfor implementation...')
    t = 0;
    tic
    AKi = cell(K,1); Bi = cell(K,1);
    for i = 1:K
        inds = (i-1)*Np + (1:Np);
        AKi{i} = AK(inds,inds);
        Bi{i} = Bh(inds,:);
    end
    t=t+toc
    tic
    parfor i = 1:K
        Tblk{i} = AKi{i}\Bi{i};
    end
    t = t+toc
end
T = cell2mat(Tblk);
A = T'*Bh;

% forcing
b = T'*M*f;

% BCs on u
u0 = zeros(size(B,2),1);

% BCs on flux
uh0 = zeros(nM,1);
bnf = nxfb*b1 + nyfb*b2; % beta_n, determines inflow vs outflow
bmaskf = xfb < -1 + NODETOL & nxfb < -NODETOL; %(bnf < NODETOL); % inflow = beta_n < 0
uh0(fmapB) = bnf.*(yfb<0).*(1+yfb);  % BC data on flux = bn*u - eps*du/dn

U0 = [u0;uh0];

% remove BCs on u on inflow for stability
vmapBTr(xr < -1+NODETOL) = [];
%   vmapBTr = []; % removes all Dirichlet BCs for testing....

% BCs on U: ordered first
b = b - A*U0;
b(vmapBTr) = U0(vmapBTr);
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
fmapB(~bmaskf) = []; % do 0 Neumann outflow BCs on test fxns

bci = nU + fmapB; % skip over u dofs
b(bci) = U0(bci);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

% Pre = buildOAS_primalDPG(Rp,A,Ntrial,Nf,fpairs,xF,yF,b1,b2);
% b = rand(nU,1);
% [U, flag, relres, iter, resvec] = fpcg(A,b,1e-6,100,@(x) Pre(x));

% b = b*0;
% b(nU + round(nM/2)) = 1;
% b(round(nU*.9)) = 1;
% keyboard
if ~plotFlag
    return
else
    U = A\b;
    u = Rr'*U(1:nU);
    plotSol(u,25);    
    title('DPG with fluxes and traces')
%     color_line3(xF,yF,U(nU+(1:nM)),U(nU+(1:nM)),'o')
end


function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
global b1
global b2
global ep
Ks = Dx'*M*Dx + Dy'*M*Dy;

S = -(b1*Dx+b2*Dy)'*M;
% S = M*(b1*Dx+b2*Dy);
Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Poisson
% Test = M + Ks;
% Trial = Ks;

% CD
Test = M + ep*Ks + Kb;
Trial = S + ep*Ks;

% Helmholtz
% k = 100;
% Test = k^2*M + Ks;
% Trial = -k^2*M + Ks;

