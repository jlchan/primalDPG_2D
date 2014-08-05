function pDPG_tensorProductOAS

%% make DPG system
Globals2D
FaceGlobals2D

% Polynomial order used for approximation
Ntrial = 10;
Ntest = 12;
Nf = 10;
% Nt = 0; % dummy

N = Ntest;
[Nv, VX, VY, K, EToV] = QuadMesh2D(1);

% Initialize solver and construct grid and metric
StartUp2D;FaceStartUp2D 

% get block operators
[M, Dx, Dy] = getBlockOps();
AK = M + Dx'*M*Dx + Dy'*M*Dy;
BK = Dx'*M*Dx + Dy'*M*Dy;

f = ones(Np*K,1);
b = M*f;

[Rp Irp vmapBTr xrB yrB xr yr] = pRestrictCG(Ntest,Ntrial); % restrict test to trial space
Rr = Rp*Irp';

Bhat = getMortarConstraint();
xfB = xf(fmapB); yfB = yf(fmapB); nxf = nxf(fmapB);nyf = nyf(fmapB);

B = BK*Rr';   % form rectangular bilinear form matrix

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
tic
for i = 1:K % independently invert
    inds = (i-1)*Np + (1:Np);
    Tblk{i} = AK(inds,inds)\Bh(inds,:);
    disp(['on element ' num2str(i)])    
end
disp(['time for test function computation = ', num2str(toc)])
T = cell2mat(Tblk);

% assemble forcing + matrix
A = T'*Bh;
b = T'*b;

% BCs on u, f
u0 = zeros(size(B,2),1);
uh0 = zeros(nM,1); 

% BCs on U: ordered first
b(vmapBTr) = 0; 
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0; A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

U = A\b;

u = Rr'*U(1:nU);
uhat = U(nU+(1:nM));

% figure
% plotFlux(uhat)
% plotSol(u,25)
% plotSol(Dx*u,25)
% plotSol(Dy*u,25) 

%% build OAS dof lists - experiment with 1D tensor productization
% [Mhandle Ak] = buildOAS_primalDPG(Rp,A,Ntrial,2); % 2 = face patch, 3 = elem patch only

[vmapC cc] = find(Rp);
vmapC = reshape(vmapC,(Ntrial+1),(Ntrial+1));

NpU = (Ntrial+1)^2;
NpF = Nfrp*4;
iU = 1:NpU;
iF = NpU+1:NpU+NpF;

[A1D M1D] = DPG_1D(Ntrial,Ntest,1);

% plot flux/field nodes
if 0
    figure
    plot(xr,yr,'.','markersize',10);hold on;
    j = 1;
    for k = 1:K
        for i = 1:size(xr,1)
            i = find(vmapC==i);
            off = .2*((k-1)/K);
            text(xr(i,k)+off,yr(i,k)+off,num2str(j),'fontsize',16)
            j = j+1;
        end        
    end
    figure
    plot(xf,yf,'x','markersize',10);hold on
    j = 1;
    for f = 1:NfacesU
        for i = 1:Nfrp
            off = .2*((f-1)/NfacesU);
            plot(xf(:,f),yf(:,f),'.-')
            text(xf(i,f)+off,yf(i,f)+off,num2str(j),'fontsize',16)
            j = j+1;
        end        
    end
end

fmap = reshape(fmap,Nfrp,Nfaces);
% u = Pre(b,vmapC,fmap,A1D,M1D); % applies 1D DPG operator in both directions
[x flag relres iter resvec] = pcg(A,b,1e-8,100,@(x) Pre(x,vmapC,fmap,A1D,M1D));
[x flag relres iter resvecNoPre] = pcg(A,b,1e-8,100);
u = Rr'*u(1:nU);
plotSol(u,25);
figure
semilogy(resvec,'.-');hold on;
semilogy(resvecNoPre,'ro-')


function b = Pre(x,vmapC,fmap,A1D,M1D)
NpTrial = size(vmapC,1);
NpU = NpTrial^2;
fmap_x = NpU + [fmap(:,4) fmap(:,2)]; % boundary mortars for one elem
fmap_y = NpU + [fmap(:,1) fmap(:,3)]; % boundary mortars for one elem
Ix = [vmapC;fmap_x'];
Iy = [vmapC fmap_y]';

% what do I do here?
b = zeros(size(x));
b(Ix) = b(Ix) + (A1D\x(Ix));
b(Iy) = b(Iy) + (A1D\x(Iy));
