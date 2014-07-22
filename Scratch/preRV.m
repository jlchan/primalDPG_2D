Globals2D
addpath('../Precond')
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell00625.neu');
% Nv = 3;
% VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
% EToV = [3 1 2];
% K = 1;

<<<<<<< HEAD
Nvec = [3 6 9];
=======
Nvec = 3:8;
>>>>>>> 6ebf5ac9b5f83ed9586d89700f7d4b84ba759ced
kvec = randperm(K);
kvec = kvec(1:3);

Nchoices = length(Nvec)*length(kvec);
C = hsv(Nchoices);
i = 1;
leg = {};
for N = Nvec; % Ntest
    
    % Initialize solver and construct grid and metric
    StartUp2D;
    
    % get block operators
    [Mh, Dxh, Dyh] = getBlockOps();
    
    % mass matrix
    MK = MassMatrix;
    GradK = [Dr;Ds]; M2K = blkdiag(MK,MK);
    AK = GradK'*blkdiag(MK,MK)*GradK;
    e = ones(Np,1)/Np;
    %     [U S V] = svd(AK);r = nnz(diag(S)>1e-8);
    %     U = U(:,1:r); D = S(1:r,1:r);
    %     P = inv(MK) - inv(MK)*U*inv(diag(1./diag(D)) + U'*inv(MK)*U)*U'*inv(MK);
    %     keyboard
    %     P = eye(size(AK));
    %     Pre = @(x) P*x;
%         Pre = @(x,h) (1/h)*MK\x + h*inv(V)'*pinv(V'*AK*V)*inv(V)*x;
    
%     Pre = @(x,h) (1/h)*MK\x + h*(1e-3*MK + GradK'*M2K*GradK)\x;
    Pre = @(x,h) (1/h)*MK\x + h*(GradK'*M2K*GradK + e*e')\x;
    for k = kvec
        inds = (k-1)*Np+(1:Np);
        M = Mh(inds,inds);Dx = Dxh(inds,inds);Dy = Dyh(inds,inds);
        Grad = [Dx;Dy];
        M2 = blkdiag(M,M);
        
        RV = M + Grad'*M2*Grad;
        
        h = J(1,k); % assumes affine
        b = rand(size(RV,1),1);
        [U, flag, relres, iter, resvec] = pcg(RV,b,1e-7,50,@(x) Pre(x,h));
        semilogy(resvec,'.-','color',C(i,:));hold on
        leg{i} = ['Matrix size ', num2str(size(RV,1))];
        i = i+1;
    end
end

legend(leg)
