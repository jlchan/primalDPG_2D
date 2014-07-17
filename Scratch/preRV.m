Globals2D
addpath('../Precond')
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell00625.neu');
% Nv = 3;
% VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
% EToV = [3 1 2];
% K = 1;

Nvec = 3:6;
kvec = randperm(K);
kvec = kvec(1:5);

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
    
    for k = kvec
        inds = (k-1)*Np+(1:Np);
        M = Mh(inds,inds);Dx = Dxh(inds,inds);Dy = Dyh(inds,inds);
        Grad = [Dx;Dy];
        M2 = blkdiag(M,M);
        
        RV = M + Grad'*M2*Grad;
        
        h = J(1,k); % assumes affine
        b = rand(size(RV,1),1);
        Pre = @(x) (1/h)*MK\x + h*(MK + GradK'*M2K*GradK)\x;
        [U, flag, relres, iter, resvec] = pcg(RV,b,1e-7,50,Pre);
        semilogy(resvec,'.-','color',C(i,:));hold on
        i = i+1;
    end
end
