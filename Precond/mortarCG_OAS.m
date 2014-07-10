function err= mortarCG_OAS(N,Nf,mesh)

Globals2D;

% Polynomial order used for approximation 
if nargin<1
    useCG = 1;
    N = 2; % when N = even, Nf = N-1, fails?
    Nf = 0;
    %     Read in Mesh 
    mesh = 'squarereg.neu';
    mesh = 'Maxwell0125.neu';
%     mesh = 'Maxwell1.neu';
% mesh = 'lshape.neu';
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
    %Nv = 3; VX = VX(EToV(1,:)); VY = VY(EToV(1,:)); EToV = [3 1 2]; K = 1;
    
else
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
end


% Initialize solver and construct grid and metric
StartUp2D;

[M, Dx, Dy] = getBlockOps();
AK = getVolOp(M,Dx,Dy);

uex = @(x,y) sin(pi*x).*sin(pi*y);
uex = project_fxn(uex,25);
f = uex*(1+2*pi^2);
f = ones(size(x(:)));

if useCG
    [R vmapBT] = getCGRestriction();        
    A = R*AK*R';
    b = R*M*f;
    
    n = size(A,1); 
    u0 = zeros(n,1); 
    % Dirichlet BC on the left
%     u0(vmapBT) = (x(vmapB) < -1+1e-7).*sqrt(1-y(vmapB).^2); 
%     vmapBT(x(vmapB) < -1+NODETOL) = []; % remove left BCs for Neumann
    b = b - A*u0;    
    b(vmapBT) = u0(vmapBT);
    A(vmapBT,:) = 0; A(:,vmapBT) = 0;
    A(vmapBT,vmapBT) = speye(length(vmapBT));
    u = A\b;    
    
    % build OAS preconditioner 
    [r c] = find(R); r = reshape(r,Np,K); c = reshape(c,Np,K);
    Ak = cell(K,1);Aki = cell(K,1);
    for k = 1:K        
%         nbr = unique([k EToE(k,:)]); % 1 elem overlap
        nbr = k; % no overlap
        inds = unique(r(:,nbr)); 
        Aki{k} = inds;
        Ak{k} = eye(numel(inds));%A(inds,inds);                   
%         clf
%         plotVerts;hold on
%         ci = unique(c(:,nbr));
%         plot(x(ci),y(ci),'ro')
%         title(['k = ', num2str(k)])
%         pause
    end    
%     u = OAS(b,Ak,Aki);
%     plotSol(R'*u,25)    
%     keyboard

    % coarse grid solver
    Rs = spdiag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
    [Rc1 Ir1 vmapBT1 xr1 yr1] = pRestrictCG(N,1); % interp down
    Rp1 = Rc1*Ir1';
    R1 = Rs*Rp1'; % interp down to P1
    
%     A1 = R1'*A*R1; 
%     A1 = Rp1*AK*Rp1';
%     A1(vmapBT1,:) = 0;A1(:,vmapBT1) = 0;
%     A1(vmapBT1,vmapBT1) = speye(numel(vmapBT1));
% %     b1 = R1'*b; b1(vmapBT1) = 0; % FIGURE OUT HOW TO USE u0!!!
%     B = speye(size(R1,2)); B(vmapBT1,vmapBT1) = 0; % imposes homogeneous BCs
%     P1 = @(b) R1*(A1\(B*R1'*b));
%     keyboard
      P1 = @(b) R1*((R1'*A*R1)\(R1'*b)); % should be able to "ignore" BC imposition - pcg acts on residual e(vmapBT) = 0
    
%     levels = agmg_setup(A1);
%     [x flag relres iter resvec1] = agmg_solve(levels, b1, 50, 1e-6);
%     semilogy(resvec1,'.-')

    [u, flag, relres, iter, resvec] = pcg(A,b,1e-6,75,@(x) OAS(x,Ak,Aki) + P1(x));   
%     [u, flag, relres, iter, resvec] = gmres(A,b,[],1e-6,50,@(x) OAS(x,Ak,Aki) + P1(x));   

semilogy(resvec,'.-'); hold on;
    title('OAS for CG')
    
    u = R'*u;    
else
%     [B, vmapBF, xfb, yfb, nxf, nyf fmap xf yf] = getMortarConstraint(Nf);    
    [B vmapBF xf yf nxf nyf fmap fpairs] = getMortarConstraint(Nf);        
    
    xfb = xf(vmapBF);yfb = yf(vmapBF);
    nxf = nxf(vmapBF);nyf = nyf(vmapBF);        
    nU = size(B,2); % num CG nodes
    nM = size(B,1); % num mortar nodes
    O = sparse(nM,nM);     
    Am = [AK B';B O];
    
    b = M*f;
    bm = [b;zeros(nM,1)];        
    um = Am\bm;
    u = um(1:Np*K);
    f = um(Np*K+1:end);
        
    % build OAS preconditioner
    S = B*(AK\B'); bS = B*(AK\b);    
    FToE = getFToE(fpairs);
    NfacesU = size(fpairs,2);
    
    fmapU = reshape(1:(Nf+1)*NfacesU,Nf+1,NfacesU);            
%         M = sparse(size(S,1),size(S,2));
    Sf = cell(NfacesU,1);
    Sfi = cell(NfacesU,1);
    for f = 1:NfacesU
%         inds = unique(fmap(:,FToE(f,:)));
        inds = fmapU(:,f); % no overlap
%                 M(inds,inds) = M(inds,inds) + inv(S(inds,inds));
        Sf{f} = S(inds,inds);
        Sfi{f} = inds;
    end    
    
    % build coarse grid solver
    If1 = ones(Nf+1,1)/sqrt(Nf+1); % interp from constant to Nf+1
    If1 = kron(speye(NfacesU),If1);    
    S1 = If1'*S*If1; bS1 = If1'*bS;
%     P1 = @(bS) If1*(S1\bS1);        
    P1 = @(bS) If1*(S1\(If1'*bS));            
%     levels1 = agmg_setup(S1);            
%     [x flag relres iter resvec1] = agmg_solve(levels1, bS1, 50, 1e-6);        
%     semilogy(resvec1,'.-')
%     keyboard
    
    [f, flag, relres, iter, resvec] = pcg(S,bS,1e-6,75,@(x) OAS(x,Sf,Sfi) + P1(x));
%     [f, flag, relres, iter, resvec] = gmres(S,bS,[],1e-6,50,@(x) OAS(x,Sf,Sfi) + P1(x));   
    semilogy(resvec,'.-');    hold on
    title(sprintf('OAS for mortars with N = %d, mesh = %s',N,mesh))    
    
    u = AK\(b-B'*f);
end

err = u-uex;
err = sqrt(err'*M*err);
return
if nargin<1
    plotSol(u,25);    
    title(sprintf('N = %d, Nf = %d, err = %d',N, Nf, err))        
end
% keyboard

function Vol = getVolOp(M,Dx,Dy)

Globals2D
Ks = Dx'*M*Dx + Dy'*M*Dy;
% b1 = 1; b2 = 1;ep = 1e-6;

% S = -(b1*Dx+b2*Dy)'*M;
% Kb = (b1*Dx+b2*Dy)'*M*(b1*Dx+b2*Dy);

% Vol = M + Kb + ep*Ks; % Convection-diffusion
% Vol = M+ 1e-4*Ks + Kb;  % Poisson
Vol = M + Ks;
% Vol = M + Kb;


function b = OAS(x,Af,Afi)
d = zeros(size(x));
for i = 1:length(Afi)
    d(Afi{i}) = d(Afi{i}) + 1;
end
d = 1./d; % scale
d = d.^0;

b = zeros(size(x));
for i = 1:length(Afi)    
    b(Afi{i}) = b(Afi{i}) +  d(Afi{i}).*(Af{i}\x(Afi{i}));    
end

