function pDPG_poisson_ref

Globals2D
FaceGlobals2D;

% Polynomial order used for approximation
Ntrial = 2;
Ntest = Ntrial + 2;
Nf = Ntrial-1;

N = Ntest;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = QuadMesh2D(16);
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('backdrop1.neu');

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

Nref = 6;
for rr = 1:Nref    
    % get block operators
    [M, Dx, Dy] = getBlockOps();
    [AK, BK] = getVolOp(M,Dx,Dy);
    f = ones(Np*K,1);
        
    [R vmapBT] = getCGRestriction();
    [Rp Irp vmapBTr xr yr] = pRestrictCG(N,Ntrial); % restrict test to trial space
    Rr = Rp*Irp';
%     Rr = Irp';
    Bhat = getMortarConstraint();
    xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);nyfb = nyf(fmapB);
    
    B = BK*Rr';  % form rectangular bilinear form matrix
    % B = BK;
    
    [nV nU] = size(B); % num test nodes, num trial nodes
    nM = size(Bhat,1); % num mortar nodes
    nTrial = nU + nM;
    
    Bh = [B Bhat'];
    Tblk = cell(K,1);    
        
    % compute trial-to-test operator
    tic
    for i = 1:K % independently invert
        inds = (i-1)*Np + (1:Np);
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
%         disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
%     disp(['time for test function computation = ', num2str(toc)])

    T = cell2mat(Tblk);
    A = T'*Bh;
        
    % forcing
    bh = M*f;
    b = T'*bh;
    
    % BCs on u
    u0 = zeros(size(B,2),1);
    
    % BCs on flux
    uh0 = zeros(nM,1);        
    U0 = [u0;uh0];
    
    % remove BCs on u on inflow for stability
    Dir = (yr > -1 + NODETOL); %| (xr > NODETOL);    
    vmapBTr(~Dir) = [];
    
    % BCs on U: ordered first
    b = b - A*U0;    
    b(vmapBTr) = U0(vmapBTr);
    A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
    A(vmapBTr,vmapBTr) = speye(length(vmapBTr));
    
    % homogeneous BCs on V implied by mortars, mortar BCs removes test
    % functions BCs.
    Neum = (abs(yfb+1) < NODETOL);% & (xfb < -NODETOL);
    fmapB(~Neum) = []; % do 0 Neumann outflow BCs on test fxns    
%     fmapB = [];
    
    bci = nU + fmapB; % skip over u dofs    
    b(bci) = uh0(fmapB);
    A(bci,:) = 0; A(:,bci)=0;
    A(bci,bci) = speye(length(bci));
        
    U = A\b;
        
    PlotMesh2D        
    if (rr < Nref) % don't refine on the last step        
        res = Bh*U-bh; 
        eK = zeros(K,1);
        eRep = zeros(Np*K,1);
        for i = 1:K % independently invert
            inds = (i-1)*Np + (1:Np);
            eRep(inds) = AK(inds,inds)\res(inds);
            eK(i) = res(inds)'*eRep(inds);            
%             disp(['on element ' num2str(i)])            
        end
        
        ref = eK > .15*max(eK);
        nRef = nnz(ref)

        %         % Dorfler
        %         theta = .2;
        %         nRef = 0;
        %         while nRef==0
        %             ref = eK > theta*sum(eK);
        %             nRef = nnz(ref);
        %             theta = theta/2;
        %         end
        
        
%         plotSol(Rr'*U(1:nU),25)
%         plotFlux(U(nU+(1:nM)));
%         plotSol(res,25)
%         plotSol(eRep,25)
%         eKplot = zeros(Np,K);
%         for k = 1:K
%             eKplot(:,k) = eK(k);
%         end
%         plotSol(eKplot,25)                
        
        % refine mesh
        eflag = zeros(K,3);
        eflag(ref,:) = 1;
        Refine2D(eflag);        
    end
end

u = Rr'*U(1:nU);
f = U(nU+(1:nM));

Nplot = 25;
plotSol(u,Nplot);
gU = (Dx*u).^2+ (Dy*u).^2;
plotSol(gU,Nplot);
plotFlux(f)
title('DPG with fluxes and traces')




function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = M + Ks;
Trial = Ks;
