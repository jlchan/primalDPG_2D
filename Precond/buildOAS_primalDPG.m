function [Mhandle Ak] = buildOAS_primalDPG(R,A,Norder,patches)

Globals2D;FaceGlobals2D

nU = size(R,1);
nM = size(A,1) - nU;

% build face data arrays
FToE = getFToE(fpairs);
NfacesU = size(fpairs,2);
Nfrp = Nf+1;
sharedFaces = ~ismember(fpairs(2,:),fpairs(1,:));

fmapU = reshape(1:(Nf+1)*NfacesU,Nf+1,NfacesU);
fmap = zeros(Nfrp,Nfaces*K);
fmap(:,fpairs(1,:)) = reshape(1:Nfrp*NfacesU,Nfrp,NfacesU);
fmap(:,fpairs(2,sharedFaces)) = fmap(:,fpairs(1,sharedFaces));
fmap = reshape(fmap,Nfrp*Nfaces,K); %fmap = elem to hybrid nodes

% build OAS preconditioner algebraically:
if (size(EToV,2)==4) % if quad
    Norderp = (Norder+1)^2;
else
   Norderp = (Norder+1)*(Norder+2)/2;
end
[rr cc] = find(R);
rr = reshape(rr,Norderp,K);cc = reshape(cc,Norderp,K);

if nargin<4
    patches = 1;
end
if patches==1
    % build patches around vertices
    vertices = rr([1 Norder+1 Norderp],:);
    vnodes = unique(vertices(:));
    Nverts = length(vnodes);
    Ak = cell(Nverts,2);
    for i = 1:length(vnodes)
        % get element patches
        [~, elems] = find(vnodes(i)==vertices);
        inds = unique(rr(:,elems));
        
        fInds = nU + fmap(:,elems); % get face inds of elems
        %     if 0 % removes inflow BCs
        %         nxf = nxf(:,elems);nyf = nyf(:,elems); % get outward normals at faces
        %         [b, ~, j] = unique(fInds(:));
        %         bNodes = accumarray(j(:),1)<=1;
        %         Binds = fInds(bNodes); % gets non-repeated elements = boundary nodes of subdomain
        %         beta_n = b1*nxf(bNodes) + b2*nyf(bNodes);
        %         inflow = beta_n < -NODETOL;
        %         Binds(inflow) = []; % inflow indices
        %         Binds = [];
        %         fInds = setdiff(fInds(:),Binds); % automatically removes non-uniquenss too...
        %     end
        fInds = unique(fInds(:));
        inds = [inds(:); fInds];
        
        %     clf;plot(x(cc(:,elems)),y(cc(:,elems)),'ms');hold on;plot(xf,yf,'.');plot(xf(inflow),yf(inflow),'ro','markersize',8)
        %     title(['min \sigma = ',num2str(min(svd(full(A(inds,inds)))))])
        %     pause
        
        Ak{i,1} = inds;
        Ak{i,2} = A(inds,inds);
        
    end
elseif patches == 2 % take face stencil    
    for i = 1:K
        % get element + face contribs
        elems = unique([i EToE(i,:)]); % face overlap only
        inds = unique(rr(:,elems));        
        fInds = nU + fmap(:,elems); % get face inds of elems
        fInds = unique(fInds(:));
        inds = [inds(:); fInds];        
        Ak{i,1} = inds;
        Ak{i,2} = A(inds,inds);
    end
else % take just a single element
    for i = 1:K
        inds = unique(rr(:,i));
        fInds = nU + fmap(:,i);
        fInds = unique(fInds(:));
        inds = [inds(:); fInds];
        Ak{i,1} = inds;
        Ak{i,2} = A(inds,inds);
    end
end

% coarse grid solver
[Rc1 Ir1 vmapBT1] = pRestrictCG(Norder,1); % interp down
Rp1 = Rc1*Ir1';
Rs = spdiag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
R1 = Rs*Rp1'; % interp down to P1
If1 = ones(Nf+1,1)/sqrt(Nf+1); % interp from constant to Nf+1
If1 = kron(speye(NfacesU),If1);

RR = blkdiag(R1,If1);
% P1 = @(b) R1*((R1'*A*R1)\(R1'*b)); % should be able to "ignore" BC imposition - pcg acts on residual e(vmapBT) = 0
% S1 = If1'*S*If1; %bS1 = If1'*bS;
% P1 = @(bS) If1*(S1\(If1'*bS));
P1 = @(b) RR*((RR'*A*RR)\(RR'*b)); % exact P1 solve

% a few AGMG sweeps instead of P1 solve
% maxiter = 5; 
% levels = agmg_setup(RR'*A*RR); 
% P1 = @(x) RR*agmg_solve(levels,RR'*x,maxiter,1e-8);

Mhandle = @(x) P1(x) + OAS(x,Ak); % coarse solver + OAS solver

