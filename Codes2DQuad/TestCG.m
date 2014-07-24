clear
% Driver script for solving the 2D vacuum Maxwell's equations on TM form
% Globals2D;
Globals2D;FaceGlobals2D;

% Polynomial order used for approximation
N = 10;
Nf = 4;
[Nv, VX, VY, K, EToV] = QuadMesh2D(1);

StartUp2D;FaceStartUp2D
useCG = 0;
[M, Dx, Dy] = getBlockOps();
Grad = [Dx;Dy];
A = M + Grad'*blkdiag(M,M)*Grad;    
f = ones(Np*K,1);

if useCG
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
    
    S = B*(A\B');
    
    [M1 Dr1 Ds1] = getFEMOps(N);
    % % scale with geom factors
    M1 = spdiag(J(:))*kron(speye(K),M1); % J = h^2
    Dx1 = spdiag(rx(:))*kron(speye(K),Dr1) + spdiag(sx(:))*kron(speye(K),Ds1);
    Dy1 = spdiag(ry(:))*kron(speye(K),Dr1) + spdiag(sy(:))*kron(speye(K),Ds1);
    
    Grad1 = [Dx1;Dy1];
    M1L = spdiag(sum(M1));
    A1 = M1 + Grad1'*blkdiag(M1L,M1L)*Grad1;
%     iA1 = diag(diag(inv(A1)));%spai(A1,.5);
    iA1 = diag(1./diag(A1));
    S1 = B*(iA1*B');
    
    spy(S);hold on;spy(S1,'ro');
    c = cond(full(S1\S));
    cs = cond(full(S));
    title(['Nnz reduction = ', num2str(nnz(S)/nnz(S1)), ', cond of S_1^{-1}S = ', num2str(c) ', vs cond of S = ', num2str(cs)])   
end
% plotSol(u,50)
