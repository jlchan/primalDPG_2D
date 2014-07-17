function Mhandle = buildOAS_mortar(S)

Globals2D;FaceGlobals2D

% build OAS preconditioner
useOverlap = 1;
if useOverlap
    FToE = getFToE(fpairs);
end
% NfacesU = size(fpairs,2); %
nM = size(S,1);
NfacesU = nM/(Nf+1);

fmapU = reshape(1:(Nf+1)*NfacesU,Nf+1,NfacesU);
%         M = sparse(size(S,1),size(S,2));
Sf = cell(NfacesU,2);
for f = 1:NfacesU
    if useOverlap
        inds = unique(fmap(:,FToE(f,:))); % 1 element overlap
    else
        inds = fmapU(:,f); % no overlap
    end
    Sf{f,1} = inds;
    Sf{f,2} = S(inds,inds);    
end

% build coarse grid solver
If1 = ones(Nfrp,1)/sqrt(Nfrp); % interp from constant to Nf+1
If1 = kron(speye(NfacesU),If1);
S1 = If1'*S*If1; %bS1 = If1'*bS;

P1 = @(bS) If1*(S1\(If1'*bS));
Mhandle = @(x) P1(x) + OAS(x,Sf); % coarse solver + OAS solver