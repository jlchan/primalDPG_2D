function primalDPG_poisson

Globals2D
FaceGlobals2D

% Polynomial order used for approximation
Ntrial = 1;
Ntest = Ntrial+1;
Nf = Ntrial-1;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
[Nv, VX, VY, K, EToV] = QuadMesh2D(2);
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D 

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = ones(Np*K,1);
b = M*f;

%[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xr yr] = pRestrictCG(Ntest,Ntrial); % restrict test to trial space
Rr = Rp*Irp';
Bhat = getMortarConstraint();
xfB = xf(fmapB); yfB = yf(fmapB); nxf = nxf(fmapB);nyf = nyf(fmapB);

B = BK*Rr';   % form rectangular bilinear form matrix

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

% penalty/robin BCs 
bmask = abs(y(vmapB)) > 1 - NODETOL; % top/bottom boundaries
[Mb Eb] = getBoundaryMatrix(bmask(:));
u0tb = 1+x(vmapB);
pen = 1;
B = B + pen*Eb'*Mb*Eb*Rr'; % this adds a penalty term on u (or Robin condition)
b = b + pen*Eb'*Mb*u0tb;

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
b = T'*b;

% BCs on u
u0 = zeros(size(B,2),1);
left = xr < -1+NODETOL;
u0(vmapBTr) = left.*sqrt(1-yr.^2); 
right = xr > 1-NODETOL;
top = yr > 1-NODETOL;
bot = yr < -1+NODETOL;
% u0(vmapBTr) = top.*sqrt(1-xr.^2);

% BCs on flux
uh0 = zeros(nM,1); 
leftf = xfB < -1+NODETOL; % right boundary
rightf = xfB > 1-NODETOL; % right boundary
uh0(fmapB) = -rightf.*nxf.*((yfB<=0) - (yfB>0));  % BC data on -du/dn
topf = yfB > 1-NODETOL;
botf = yfB < -1+NODETOL;
U0 = [u0;uh0];

b = b - A*U0; % get lift

% BCs on U: ordered first
vmapBTr(~left) = [];
b(vmapBTr) = U0(vmapBTr);
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
fmapB(leftf) = [];
% fmapB = [];
bci = nU + fmapB; % skip over u dofs
b(bci) = uh0(fmapB);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

U = A\b;
u = Rr'*U(1:nU);
% uhat = U(nU+(1:nM));
% figure
% color_line3(xf,yf,uhat,uhat,'.');

% Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot); 
plotSol(u,55)
keyboard
% color_line3(x,y,u,u,'.')
title('DPG with fluxes and traces')


function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = M + Ks;
Trial = Ks;

