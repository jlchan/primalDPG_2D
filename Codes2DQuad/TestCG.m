clear
% Driver script for solving the 2D vacuum Maxwell's equations on TM form
% Globals2D;
Globals2D;FaceGlobals2D;

% Polynomial order used for approximation
N = 4;
Nf = 0;
[Nv, VX, VY, K, EToV] = QuadMesh2D(4);

StartUp2D;FaceStartUp2D
CG = 0;
[M, Dx, Dy] = getBlockOps();
Grad = [Dx;Dy];
A = M + Grad'*blkdiag(M,M)*Grad;    
f = ones(Np*K,1);

if CG
    [R vmapBT] = getCGRestriction();
    
    A = R*A*R';
    b = R*M*f;
    
    n = size(A,1);
    b(vmapBT) = 0;
    A(vmapBT,:) = 0; A(:,vmapBT) = 0;
    A(vmapBT,vmapBT) = speye(length(vmapBT));
    u = A\b;
    u = R'*u;
else
    B = getMortarConstraint();
    
    nU = Np*K; % num CG nodes
    nM = Nfrp*NfacesU; % num mortar nodes
    O = sparse(nM,nM);
    Am = [A B';B O];
    
    b = M*f;
    bm = [b;zeros(nM,1)];
    
    um = Am\bm;
    u = um(1:Np*K);
    f = um(Np*K+1:end);
end
plotSol(u,50)
