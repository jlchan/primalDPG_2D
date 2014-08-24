function primalDPG_poisson

Globals2D
FaceGlobals2D

% Polynomial order used for approximation
Ntrial = 2;
Ntest = Ntrial+1;
Nf = Ntrial-1;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');

[Nv, VX, VY, K, EToV] = QuadMesh2D(4);

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D 

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = ones(Np*K,1);
b = M*f;

NpTest = Np;

%[R vmapBT] = getCGRestriction();
Bhat = getMortarConstraint();

% penalty/robin BCs 
bmask = abs(y(vmapB)) > 1 - NODETOL; % top/bottom boundaries
[Mb Eb] = getBoundaryMatrix(bmask(:));
u0tb = 1+x(vmapB);
pen = 1;
BK = BK + pen*Eb'*Mb*Eb; % this adds a penalty term on u (or Robin condition)
b = b + pen*Eb'*Mb*u0tb;

% switches to Ntrial globals
[Rp Irp vmapBT] = pRestrictCG(Ntest,Ntrial); % restrict test to trial space
Rr = Rp*Irp';
B = BK*Rr';   % form rectangular bilinear form matrix
[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
tic
for i = 1:K % independently invert
    inds = (i-1)*NpTest + (1:NpTest);
    Tblk{i} = AK(inds,inds)\Bh(inds,:);
    disp(['on element ' num2str(i)])
    %             Tblk{i} = Bh(inds,:);
end
disp(['time for test function computation = ', num2str(toc)])
T = cell2mat(Tblk);
A = T'*Bh;

% forcing
b = T'*b;

% BCs on u
u0 = zeros(size(B,2),1);
left = x(vmapB) < -1+NODETOL;
u0(vmapBT) = left.*sqrt(1-y(vmapB).^2); 

% BCs on flux
xfB = xf(fmapB); yfB = yf(fmapB); nxf = nxf(fmapB);nyf = nyf(fmapB);
uh0 = zeros(nM,1); 
leftf = xfB < -1+NODETOL; % right boundary
rightf = xfB > 1-NODETOL; % right boundary
uh0(fmapB) = -rightf.*nxf.*((yfB<=0) - (yfB>0));  % BC data on -du/dn
U0 = [u0;uh0];

b = b - A*U0; % get lift

% BCs on U: ordered first
vmapBT(~left) = [];
b(vmapBT) = U0(vmapBT);
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
fmapB(leftf) = [];
bci = nU + fmapB; % skip over u dofs
b(bci) = uh0(fmapB);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

U = A\b;
u = Rp'*U(1:nU);

plotSol(u,25)

title('DPG with fluxes and traces')


function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = M + Ks;
Trial = Ks;


