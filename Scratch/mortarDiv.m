function mortarDiv

Globals2D;

N = 6; % when N = even, Nf = N-1, fails?
Nf = 4;
%     Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
%     Nv = 3;
%     VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
%     EToV = [3 1 2];
%     K = 1;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
%     [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');

% Initialize solver and construct grid and metric
StartUp2D;

[M, Dx, Dy] = getBlockOps();
Div = [Dx Dy];
M2 = blkdiag(M,M);
A = M2 + Div'*M*Div;
f = ones(size(x(:)));
f = zeros(size(x)); f(:,K-1) = 1; f= f(:);
f2 = M2*[f;0*f];

[Bdiv vmapBF xf yf nxf nyf fpairs] = getMortarConstraintDiv(Nf);
 
xfb = xf(vmapBF);yfb = yf(vmapBF);
nxf = nxf(vmapBF);nyf = nyf(vmapBF);
nU = size(Bdiv,2); % num CG nodes
nM = size(Bdiv,1); % num mortar nodes
O = sparse(nM,nM);
Am = [A Bdiv';Bdiv O];
b2 = M2*[f;f];
bm = [b2;zeros(nM,1)];
% keyboard
um = Am\bm;

U = um(1:nU);

Iu1 = 1:Np*K;
Iu2 = Np*K + (1:Np*K);

u1 = U(Iu1);
u2 = U(Iu2);
unorm = sqrt(u1.^2 + u2.^2);
plotSol(unorm,25);
plotSol(u1,25);
plotSol(u2,25);
