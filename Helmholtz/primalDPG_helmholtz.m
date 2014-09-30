function primalDPG_helmholtz

Globals2D
FaceGlobals2D

% Polynomial order used for approximation
Ntrial = 3;
Ntest = Ntrial+1;
Nf = Ntrial-1;

N = Ntest;

global k;
k = 25; x0 = [-2.5, -.25];

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');

% [Nv, VX, VY, K, EToV] = QuadMesh2D(8);

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D

Nref = 5;
for rr = 1:Nref
    
    % get block operators
    [M, Dx, Dy] = getBlockOps();
    [AK, BK] = getVolOp(M,Dx,Dy);
    f = zeros(Np*K,1);
    b = M*f;
    
    NpTest = Np;
    
    %[R vmapBT] = getCGRestriction();
    Bhat = getMortarConstraint();
    
    % penalty/robin BCs
    bmask = ones(size(vmapB));%abs(y(vmapB)) > 1 - NODETOL; % top/bottom boundaries
    [Mb Eb] = getBoundaryMatrix(bmask(:));
    BK = BK + 1i*k*Eb'*Mb*Eb; % this adds a penalty term on u (or Robin condition)
    % BK = BK + k*Eb'*Mb*Eb; % this adds a penalty term on u
    % switches to Ntrial globals
    
    [Rp Irp vmapBT] = pRestrictCG(Ntest,Ntrial); % restrict test to trial space
    Rr = Rp*Irp';

    B = BK*Rr';   % form rectangular bilinear form matrix
    [nV nU] = size(B); % num test nodes, num trial nodes
    nM = size(Bhat,1); % num mortar nodes
    nTrial = nU + nM;
    
    Bh = [B Bhat'];
    Tblk = cell(K,1);
    tic
    for i = 1:K % independently invert
        inds = (i-1)*NpTest + (1:NpTest);
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
        disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
    disp(['time for test function computation = ', num2str(toc)])
    T = cell2mat(Tblk);
    A = T'*Bh;
    
    % forcing
    b = T'*b;
    
    % BCs on u
    u0 = zeros(size(B,2),1);
%     u0(vmapBT) = exact_sol(x(vmapB),y(vmapB),k,x0);
    
    % BCs on flux
    xfB = xf(fmapB); yfB = yf(fmapB); nxf = nxf(fmapB);nyf = nyf(fmapB);
    uh0 = zeros(nM,1);
    
    [u Dxu Dyu] = exact_sol(xfB,yfB,k,x0);
    uh0(fmapB) = nxf.*Dxu + nyf.*Dyu - 1i*k*u;
    
    U0 = [u0;uh0];
        
    b = b - A*U0; % get lift
    
    % BCs on U: ordered first
%     b(vmapBT) = U0(vmapBT);
%     A(vmapBT,:) = 0; A(:,vmapBT) = 0;
%     A(vmapBT,vmapBT) = speye(length(vmapBT));
    
    % homogeneous BCs on V are implied by mortars.
    % BCs on mortars removes BCs on test functions.
    % fmapB(leftf) = [];
    bci = nU + fmapB; % skip over u dofs
    b(bci) = uh0(fmapB);
    A(bci,:) = 0; A(:,bci)=0;
    A(bci,bci) = speye(length(bci));
    
    U = A\b;
    
    if (rr < Nref) % don't refine on the last step
        
        res = Bh*U-Bh*U0;
        eK = zeros(K,1);
        eRep = zeros(Np*K,1);
        for i = 1:K % independently invert
            inds = (i-1)*NpTest + (1:NpTest);
            eRep(inds) = AK(inds,inds)\res(inds);
            eK(i) = abs(res(inds)'*eRep(inds));
            disp(['error on element ' num2str(i) ' = ', num2str(eK(i))])
        end
        ref = eK > .2*max(eK);
%         
        PlotMesh2D
%         
        % refine mesh
        eflag = zeros(K,3);
        eflag(ref,:) = 1;
               
        
        N = Ntest;
        StartUp2D;FaceStartUp2D;
%         
        Refine2D(eflag);
        
    end
end

u = Rp'*U(1:nU);

[cubR,cubS,cubW,Ncub] = Cubature2D(k);
[xcub ycub] = getGlobalNodes(cubR,cubS);
Vcub = Vandermonde2D(N,cubR,cubS);
Interp = Vcub*invV; % interp to cubature points

uexact = exact_sol(xcub,ycub,k,x0);
up = (V*V')*Interp'*diag(cubW)*reshape(uexact,Ncub,K); up = up(:);

err = reshape(uexact,Ncub,K)-Interp*reshape(u,Np,K);
L2err = diag(sqrt(cubW))*err;
L2err = norm(L2err(:));

errP = reshape(uexact,Ncub,K)-Interp*reshape(up,Np,K);
L2errP = diag(sqrt(cubW))*errP;
L2errP = norm(L2errP(:));

plotSol(u,50)
% plotSol(u-exact_sol(x,y,k,x0),25)
title(['Error = ' num2str(L2err), ', L2 proj err = ', num2str(L2errP)])
%title('DPG with fluxes and traces')

return

[M Ak] = buildOAS_primalDPG(Rp,A,Ntrial);
[U, flag, relres, iter, resvec] = gmres(A,b,100,1e-6,100,@(x) M(x));
name = ['k = ' num2str(k) ', p = ', num2str(Ntrial), ', ', num2str(K), ' elements'];
figure(2)
NtrialMax = 10;
C = hsv(NtrialMax);
semilogy(resvec,'.-','color',C(Ntrial,:),'DisplayName',name);hold on
xlabel('Number of iterations')
legend('off');legend('show')

% keyboard

function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
global k
Ks = Dx'*M*Dx + Dy'*M*Dy;

% Poisson
Test = k^2*M + Ks;
Trial = k^2*M-Ks;

function [u Dxu Dyu] = exact_sol(x,y,k,x0)

% Hankel
r = sqrt((x-x0(1)).^2+(y-x0(2)).^2);
drx = (x-x0(1))./r;
dry = (y-x0(2))./r;
u = besselh(0,k*r(:));
Dxu = .5*k*drx.*(besselh(-1,k*r) - besselh(1,k*r));
Dyu = .5*k*dry.*(besselh(-1,k*r) - besselh(1,k*r));

% Beam

