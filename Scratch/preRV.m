Globals2D
addpath('../Precond')
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell00625.neu');
% Nv = 3;
% VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
% EToV = [3 1 2];
% K = 1;

Nvec = [3 4 10 16];
kvec = randperm(K);
kvec = kvec(1:3);

Nchoices = length(Nvec)*length(kvec);
C = hsv(Nchoices);
i = 1;
for N = Nvec; % Ntest
    
    % Initialize solver and construct grid and metric
    StartUp2D;
    
    % get block operators
    [Mh, Dxh, Dyh] = getBlockOps();
    
    % mass matrix
    MK = MassMatrix;
    GradK = [Dr;Ds]; M2K = blkdiag(MK,MK);
    AK = GradK'*M2K*GradK;
    %     e = ones(Np,1)/Np;
    [U S V] = svd(AK);r = nnz(diag(S)>1e-8);
    Z = U(:,r+1:end);
    Pre = @(x,h) (1/h)*MK\x + h*(AK + Z*Z')\x;

            
    for k = kvec
        inds = (k-1)*Np+(1:Np);
        M = Mh(inds,inds);Dx = Dxh(inds,inds);Dy = Dyh(inds,inds);
                
        RVK = MK + AK;
        iRV = spai(RVK,.1);
        
        [M1 Dx1 Dy1] = getFEMOps(N);%x(:,k),y(:,k));
        Grad1 = [Dx1;Dy1];
        RV1K = M1 + Grad1'*blkdiag(M1,M1)*Grad1;               
        iRV1 = spai(RV1K,.01);
        
        return
        
        Grad = [Dx;Dy];
        M2 = blkdiag(M,M);
        
        RV = M + Grad'*M2*Grad;
        
        h = J(1,k); % assumes affine
        b = rand(size(RV,1),1);
        [U, flag, relres, iter, resvec] = pcg(RV,b,1e-7,50,@(x) Pre(x,h));
        semilogy(resvec,'.-','color',C(i,:));hold on
        i = i+1;
    end
end
