function [Mhandle Ak] = buildOAS_primalDPG(R,A,Norder)

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
    fInds = unique(fInds(:));
    inds = [inds(:); fInds];
        
    Ak{i,1} = inds;
    Ak{i,2} = A(inds,inds);
    
end

% coarse grid solver
[Rc1 Ir1 vmapBT1] = pRestrictCG(Norder,1,true); % interp down
Rp1 = Rc1*Ir1';
Rs = spdiag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
R1 = Rs*Rp1'; % interp down to P1
If1 = ones(Nf+1,1)/sqrt(Nf+1); % interp from constant to Nf+1
If1 = kron(speye(NfacesU),If1);

RR = blkdiag(R1,If1);
P1 = @(b) R1*((R1'*A*R1)\(R1'*b)); % should be able to "ignore" BC imposition - pcg acts on residual e(vmapBT) = 0
% S1 = If1'*S*If1; %bS1 = If1'*bS;
% P1 = @(bS) If1*(S1\(If1'*bS));
P1 = @(b) RR*((RR'*A*RR)\(RR'*b));

Mhandle = @(x) P1(x) + OAS(x,Ak); % coarse solver + OAS solver

function b = OAS(x,Af)
nSubD = size(Af,1);

b = zeros(size(x));
for i = 1:nSubD
    b(Af{i,1}) = b(Af{i,1}) +  Af{i,2}\x(Af{i,1}); % shouldn't scale overlapping contributions
end
