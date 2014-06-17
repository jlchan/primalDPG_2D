function [A b nU nM Np bK Rp Irp, M] = primalDPG_poisson(mesh,Ntrial,Ntest,Nflux,plotFlag)

Globals2D

N = Ntest;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);

% Initialize solver and construct grid and metric
StartUp2D;

% get block operators
[M, Dx, Dy] = getBlockOps();
[AK, BK] = getVolOp(M,Dx,Dy);
f = ones(Np*K,1);
b = M*f;

[R vmapBT] = getCGRestriction();
[Rp Irp vmapBTr xr yr] = pRestrictCG(Ntrial); % restrict test to trial space
Rr = Rp*Irp';
[Bhat vmapBF xf yf nxf nyf] = getMortarConstraint(Nflux);

B = BK*Rr';   % form rectangular bilinear form matrix

[nV nU] = size(B); % num test nodes, num trial nodes
nM = size(Bhat,1); % num mortar nodes
nTrial = nU + nM;

% penalty/robin BCs 
bmask = abs(y(vmapB)) > 1 - NODETOL; % top/bottom boundaries
[Mb Eb] = getBoundaryMatrix(bmask(:));
u0tb = 1+x(vmapB);
% B = B + 1e6*Eb'*Mb*Eb*Rr'; % this adds a penalty term on u (or Robin condition)
% b = b + 1e6*Eb'*Mb*u0tb;

Bh = [B Bhat'];
Tblk = cell(K,1);
bK = cell(K,1);

% interpolate down to 
Ir = InterpDown(Ntrial);

if 1
    tic
    for i = 1:K % independently invert
        inds = (i-1)*Np + (1:Np);        
        Tblk{i} = AK(inds,inds)\Bh(inds,:);
        disp(['on element ' num2str(i)])
        bKr = BK(inds,inds)*Ir;
        bK{i} = bKr'*(AK(inds,inds)\bKr);        
    end
    bK = blkdiag(bK{:});
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
    disp(['time for serial storage = ', num2str(toc)])    
    myCluster = parpool('local');
    tic
    parfor i = 1:K
        Tblk{i} = AKi{i}\Bi{i};        
    end
    disp(['time for test function computation = ', num2str(toc)])
end
T = cell2mat(Tblk);
A = T'*Bh;

% forcing
b = T'*b;

% BCs on u
u0 = zeros(size(B,2),1);
left = xr < -1+NODETOL;
u0(vmapBTr) = left.*sqrt(1-yr.^2); 
right = xr > 1-NODETOL;
top = yr > 1-NODETOL;
bot = yr < -1+NODETOL;
% u0(vmapBTr) = top.*sqrt(1-xr.^2);

% BCs on flux
uh0 = zeros(nM,1); 
leftf = xf < -1+NODETOL; % right boundary
rightf = xf > 1-NODETOL; % right boundary
uh0(vmapBF) = -rightf.*nxf.*((yf<=0) - (yf>0));  % BC data on -du/dn
topf = yf > 1-NODETOL;
botf = yf < -1+NODETOL;
U0 = [u0;uh0];

b = b - A*U0; % get lift

% BCs on U: ordered first
vmapBTr(~left) = [];
b(vmapBTr) = U0(vmapBTr);
A(vmapBTr,:) = 0; A(:,vmapBTr) = 0;
A(vmapBTr,vmapBTr) = speye(length(vmapBTr));

% homogeneous BCs on V are implied by mortars.
% BCs on mortars removes BCs on test functions.
vmapBF(leftf) = [];
% vmapBF = [];
bci = nU + vmapBF; % skip over u dofs
b(bci) = uh0(vmapBF);
A(bci,:) = 0; A(:,bci)=0;
A(bci,bci) = speye(length(bci));

if ~plotFlag
    return 
else
    U = A\b;
    u = Rr'*U(1:nU);
    
    Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot);
    % Nplot = 25; [xu,yu] = EquiNodes2D(Nplot);
    [ru, su] = xytors(xu,yu);
    Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
    xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
    yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
    figure
    color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
    
    title('DPG with fluxes and traces')
end

% keyboard

function [Test, Trial] = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;
Kb = (Dx+Dy)'*M*(Dx+Dy);
S = -(Dx+Dy)'*M;

% Poisson
Test = M + Ks;
Trial = .01*M+Ks;

% CD
% ep = .1;
% Test = M + ep*Ks + Kb;
% Trial = ep*Ks + S;

% % Helmholtz
% k = 50;
% Test = k^2*M + Ks;
% Trial = k^2*M-Ks;
% Test = Trial*Trial' + 1e-7*M;%speye(size(M));


function [M, Dx, Dy] = getBlockOps()

Globals2D

blkDr = kron(speye(K),Dr);
blkDs = kron(speye(K),Ds);
blkM = kron(speye(K),MassMatrix);

M = spdiag(J(:))*blkM; % J = h^2
Dx = spdiag(rx(:))*blkDr + spdiag(sx(:))*blkDs;
Dy = spdiag(ry(:))*blkDr + spdiag(sy(:))*blkDs;


function Ir = InterpDown(Nr)

Globals2D
[xr,yr] = Nodes2D(Nr); [rr, sr] = xytors(xr,yr);
Vr = Vandermonde2D(Nr,rr,sr); Ir = Vandermonde2D(Nr,r,s)/Vr;  % interp from test to trial
