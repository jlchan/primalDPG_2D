function CD_splitting
Nvec = 2:5;
mesh = {'Maxwell05.neu','Maxwell025.neu'}%,'Maxwell0125.neu'}%,'Maxwell00625.neu'};
C = hsv(length(mesh)*length(Nvec));

k = 1;
leg = {};
for i = 1:length(Nvec)
    N = Nvec(i);
    for m = 1:length(mesh)
        run(N,mesh{m},C(k,:));
        leg{k} = ['N = ' num2str(N) ' on mesh ' num2str(m)];
        k = k+1;
    end
end
legend(leg);


function run(Ntrial,mesh,c)

Globals2D

% Polynomial order used for approximation
% Ntrial = 6;
Ntest = Ntrial + 2;
% Nflux = Ntrial;
Nflux = 0;

N = Ntest;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('lshape.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell025.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell0125.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);

% Initialize solver and construct grid and metric
StartUp2D;

global b1
global b2
global ep
b1 = 1; b2 = 0;ep = 0e-2;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = 0*ones(Np*K,1);
% f = y(:)<=0;
% f = sin(pi*x(:)).*sin(pi*y(:));

[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xr yr] = pRestrictCG(Ntrial); % restrict test to trial space
Rr = Rp*Irp';
Rr = Irp';
[Bhat vmapBF xf yf nxf nyf] = getMortarConstraint(Nflux);

B = BK*Rr';   % form rectangular bilinear form matrix
% B = BK;

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

Bh = [B Bhat'];
Tblk = cell(K,1);
if 1
    tic
    for i = 1:K % independently invert
        inds = (i-1)*Np + (1:Np);
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
        disp(['on element ' num2str(i)])
        %             Tblk{i} = Bh(inds,:);
    end
    disp(['time for test function computation = ', num2str(toc)])
else
    disp('parfor implementation...')
    t = 0;
    tic
    AKi = cell(K,1); Bi = cell(K,1);
    for i = 1:K
        inds = (i-1)*Np + (1:Np);
        AKi{i} = AK(inds,inds);
        Bi{i} = Bh(inds,:);
    end
    t=t+toc
    tic
    parfor i = 1:K
        Tblk{i} = AKi{i}\Bi{i};
    end
    t = t+toc
end
T = cell2mat(Tblk);
A = T'*Bh;

% forcing
b = T'*M*f;

% BCs on u
u0 = zeros(size(B,2),1);

% BCs on flux
uh0 = zeros(nM,1);
bnf = nxf*b1 + nyf*b2; % beta_n, determines inflow vs outflow
bmaskf = (bnf < NODETOL); % inflow = beta_n < 0
uh0(vmapBF) = bnf.*(yf<0).*(1+yf);  % BC data on flux = bn*u - eps*du/dn

U0 = [u0;uh0];

% remove BCs on u on inflow for stability
% vmapBTr(xr < -1+NODETOL) = [];
  vmapBTr = []; % removes all Dirichlet BCs for testing....

% BCs on U: ordered first
b = b - A*U0;
b(vmapBTr) = U0(vmapBTr);
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
vmapBF(bmaskf) = []; % do 0 Neumann outflow BCs on test fxns

bci = nU + vmapBF; % skip over u dofs
b(bci) = uh0(vmapBF);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

useDD = 1;
if useDD
    Uex = A\b;
    uex = Uex(1:nU);
    uh = zeros(nM,1);
    nIter = 50;
    r = zeros(nIter,1);
    for k = 1:nIter
        I1 = 1:nU;
        I2 = nU + (1:nM);
        A1 = A(I1,I1);
        B = A(I1,I2);
        C = A(I2,I2);
        b1 = b(I1);
        b2 = b(I2);
%         Mr = Rr*M*Rr';        
%         D = blkdiag(A1,C);

        u =  A1\(b1-B*uh);
        uh = C\(b2-B'*u); 
%         levelsA1 = agmg_setup(A1);
%         levelsC = agmg_setup(C);
%         [u flag relres iter resvec] = agmg_solve(levelsA1, b1-B*uh, 10*N, 1e-6);
%         [uh flag relres iter resvec] = agmg_solve(levelsC, b2-B'*u, 250, 1e-8);        
        r(k) = norm(u-uex); 
        
%         u = Rr'*u;
%         Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot);
%         Nplot = 25; [xu,yu] = EquiNodes2D(Nplot);
%         [ru, su] = xytors(xu,yu);
%         Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
%         xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
%         yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
%         color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
%         title(['k = ' num2str(k)])
%         pause        
    end          
    semilogy(r,'.-','color',c);hold on;
%     keyboard
else
    U = A\b;
    u = Rr'*U(1:nU);

    
    %     color_line3(x,y,u,u,'.');
    %     return
    
    Nplot = 25;
    [xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
    Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
    xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
    yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
    figure
    color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
    
    title('DPG with fluxes and traces')
    keyboard
end


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

