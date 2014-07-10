function Mhandle = buildOAS_mortar(S,Nf,fpairs)

Globals2D;

% build OAS preconditioner
useOverlap = 0;
if nargin==3
    FToE = getFToE(fpairs);
    useOverlap = 1;
else
    FToE = [];
end
% NfacesU = size(fpairs,2); %
nM = size(S,1);
NfacesU = nM/(Nf+1);

fmapU = reshape(1:(Nf+1)*NfacesU,Nf+1,NfacesU);
%         M = sparse(size(S,1),size(S,2));
Sf = cell(NfacesU,1);
Sfi = cell(NfacesU,1);
for f = 1:NfacesU
    if useOverlap
        inds = unique(fmap(:,FToE(f,:))); % 1 element overlap
    else
        inds = fmapU(:,f); % no overlap
    end
    Sf{f} = S(inds,inds);
    Sfi{f} = inds;
end

% build coarse grid solver
If1 = ones(Nf+1,1)/sqrt(Nf+1); % interp from constant to Nf+1
If1 = kron(speye(NfacesU),If1);
S1 = If1'*S*If1; %bS1 = If1'*bS;

P1 = @(bS) If1*(S1\(If1'*bS));
Mhandle = @(x) P1(x) + OAS(x,Sf,Sfi); % coarse solver + OAS solver