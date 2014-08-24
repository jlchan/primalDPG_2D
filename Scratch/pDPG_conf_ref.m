function pDPG_conf_ref

Globals2D
FaceGlobals2D;

% Polynomial order used for approximation
Ntrial = 2;
Ntest = Ntrial + 2;
Nf = Ntrial;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
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

global b1
global b2
global ep
b1 = 1; b2 = 0;ep = 1e-4;

Nref = 7;
for rr = 1:Nref    
    % get block operators
    [M, Dx, Dy] = getBlockOps();
    [AK, BK] = getVolOp(M,Dx,Dy);
%     f = ones(Np*K,1);
    
    f = y(:)<=0;
%     f = sin(pi*x(:)).*sin(pi*y(:));
    
    [R vmapBT] = getCGRestriction();
    [Rp Irp vmapBTr xr yr] = pRestrictCG(N,Ntrial); % restrict test to trial space
    Rr = Rp*Irp';
%     Rr = Irp';
    Bhat = getMortarConstraint();
    xfb = xf(fmapB); yfb = yf(fmapB); nxfb = nxf(fmapB);nyfb = nyf(fmapB);
    
    B = BK*Rr';   % form rectangular bilinear form matrix
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
        disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
    disp(['time for test function computation = ', num2str(toc)])

    T = cell2mat(Tblk);
    A = T'*Bh;
        
    % forcing
    bh = M*f;
    b = T'*bh;
    
    % BCs on u
    u0 = zeros(size(B,2),1);
    
    % BCs on flux
    uh0 = zeros(nM,1);
    bnf = nxfb*b1 + nyfb*b2; % beta_n, determines inflow vs outflow
    inflow = (bnf < -NODETOL); % inflow = beta_n < 0
    uh0(fmapB) = 0*inflow.*bnf.*(yfb<0).*(1+yfb);  % BC data on flux = bn*u - eps*du/dn    
%     uh0(fmapB) = inflow.*bnf.*sin(yfb*pi*2);  % BC data on flux = bn*u - eps*du/dn    
    U0 = [u0;uh0];
    
    % remove BCs on u on inflow for stability
    vmapBTr(xr < -1+NODETOL) = [];
    % wall = (abs(yr+1)<NODETOL) & (xr > -NODETOL);
    % vmapBTr(~wall) = [];
%     vmapBTr = []; % removes all Dirichlet BCs for testing....
    
    % BCs on U: ordered first
    b = b - A*U0;    
    b(vmapBTr) = U0(vmapBTr);
    A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
    A(vmapBTr,vmapBTr) = speye(length(vmapBTr));
    
    % homogeneous BCs on V implied by mortars, mortar BCs removes test
    % functions BCs.
    fmapB(~inflow) = []; % do 0 Neumann outflow BCs on test fxns    
%     fmapB = [];
    
    bci = nU + fmapB; % skip over u dofs    
    b(bci) = uh0(fmapB);
    A(bci,:) = 0; A(:,bci)=0;
    A(bci,bci) = speye(length(bci));
        
    U = A\b;
    
    if (rr < Nref) % don't refine on the last step
        
        res = Bh*U-bh;
        eK = zeros(K,1);
        eRep = zeros(Np*K,1);
        for i = 1:K % independently invert
            inds = (i-1)*Np + (1:Np);
            eRep(inds) = AK(inds,inds)\res(inds);
            eK(i) = res(inds)'*eRep(inds);            
            disp(['on element ' num2str(i)])            
        end        
        ref = eK > .25*max(eK);
        
%         plotSol(Rr'*U(1:nU),25)
%         plotFlux(U(nU+(1:nM)));
%         plotSol(res,25)
%         plotSol(eRep,25)
%         eKplot = zeros(Np,K);
%         for k = 1:K
%             eKplot(:,k) = eK(k);
%         end
%         plotSol(eKplot,25)
        
        % additional regularization based on error rep
        if ((rr > 1) && false)
            Grad = [Dx;Dy];            
            eKD = kron(spdiag(eK),spdiag(ones(Np,1)));            
            reg = eKD*Grad'*blkdiag(M,M)*Grad;            
            A(1:nU,1:nU) = A(1:nU,1:nU) + 10*Rr*reg*Rr';
            plotSol(Rr'*U(1:nU),25);title('old sol')
            
            U = A\b;            
            plotSol(Rr'*U(1:nU),25);title('with regularization')
        end
        PlotMesh2D        
        
        % refine mesh
        eflag = zeros(K,3);
        eflag(ref,:) = 1;
        Refine2D(eflag);
        
    end
end

u = Rr'*U(1:nU);
f = U(nU+(1:nM));
%     color_line3(x,y,u,u,'.');
%     return

Nplot = 25;
plotSol(u,Nplot);
plotFlux(f)
title('DPG with fluxes and traces')




function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
global b1
global b2
global ep
Ks = Dx'*M*Dx + Dy'*M*Dy;

S = -(b1*Dx+b2*Dy)'*M;
% S = M*(b1*Dx+b2*Dy);
Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Poisson
% Test = M + Ks;
% Trial = Ks;

% CD
Test = M + ep*Ks + Kb;
Trial = S + ep*Ks;

% Helmholtz
% k = 100;
% Test = k^2*M + Ks;
% Trial = -k^2*M + Ks;

