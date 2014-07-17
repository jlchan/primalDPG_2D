Globals2D;

% Polynomial order used for approximation
N = 5; % when N = even, Nf = N-1, fails?
Nf = 3;
%     Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
%     Nv = 3; VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
%     EToV = [3 1 2]; K = 1;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');

% Initialize solver and construct grid and metric
StartUp2D;
FaceStartUp2D;

[M, Dx, Dy] = getBlockOps();
Grad = [Dx;Dy];M2 = blkdiag(M,M);
A = M + Grad'*M2*Grad;

% map field to interface dofs
NfpT = numel(fM); % or fP...
Efm = sparse(1:NfpT,fM(:),1,NfpT,Np*K);
Efp = sparse(1:NfpT,fP(:),fM~=fP,NfpT,Np*K); % don't subtract off vmapP for boundary nodes, where fM==fP
Ef = Efm-Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms

% make edge integration matrix on non-boundary edges
r1D = JacobiGL(0,0,N); V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

% get flux ids on boundaries
if (Nf==0)
    R1D = ones(Nfp,1)/sqrt(Nfp); % interpolate constant to multiple nodes
else
    V1Dfr = Vandermonde1D(Nf,rfr);
    V1Df = Vandermonde1D(Nf,r1D);
    R1D = V1Df/V1Dfr; % interpolate lower order polynomial to order Nf polynomial
end
sJReduc = reshape(sJ,Nfp,Nfaces*K);
sJReduc = sJReduc(:,fpairs(1,:)); % get unique faces
SJ = spdiag(sJReduc(:));
B = kron(speye(NfacesU),R1D'*M1D)*SJ*Ef; % for fluxes

f = ones(size(x(:)));

nU = Np*K; % num CG nodes
nM = NfacesU*Nfrp; % num mortar nodes
O = sparse(nM,nM);
Am = [A B';B O];

b = M*f;
bm = [b;zeros(nM,1)];

um = Am\bm;
u = um(1:Np*K);
f = um(Np*K+1:end);

plotSol(u,25);
% plotVerts
hold on
color_line3(xf,yf,f,f,'o')


