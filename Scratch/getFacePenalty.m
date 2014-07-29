
function [Bpen] = getFacePenalty()


Globals2D

fM = reshape(vmapM,Nfp,Nfaces*K); fP = reshape(vmapP,Nfp,Nfaces*K); 

% sort by node id in face and then get unique pairs of columns 
[tf, loc] = ismember(sort(fM,1)',sort(fP,1)','rows'); 
fpairs = [(1:length(loc))' loc(:)];
fpairs = unique(sort(fpairs,2),'rows')'; 
NfacesU = size(fpairs,2);

% fpairs(2,:) are duplicated face nodes, ignore them
fM = fM(:,fpairs(1,:)); fP = fP(:,fpairs(1,:));

% map field to interface dofs
NfpT = numel(fM); % or fP...
Efm = sparse(1:NfpT,fM(:),fM~=fP,NfpT,Np*K);
Efp = sparse(1:NfpT,fP(:),fM~=fP,NfpT,Np*K); % don't subtract off vmapP for boundary nodes, where fM==fP
Efm = Efm-Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms

% make edge integration matrix on non-boundary edges
r1D = JacobiGL(0,0,N); V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

sJReduc = reshape(sJ,Nfp,Nfaces*K); 
sJReduc = sJReduc(:,fpairs(1,:)); % get unique faces
SJ = spdiag(sJReduc(:));
Bpen = Efm'*kron(speye(NfacesU),M1D)*SJ*Efm; % for fluxes
