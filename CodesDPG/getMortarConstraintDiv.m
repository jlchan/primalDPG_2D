% returns constraint matrix for Hdiv vector elements in (qx,qy) ordering
% currently using discontinuous traces (should be continuous)

function B = getMortarConstraintDiv()

Globals2D;FaceGlobals2D;

[fM fP fpairs] = getFaceInfo();

% map field to interface dofs
NfpT = numel(fM); % or fP...
Efm = sparse(1:NfpT,fM(:),1,NfpT,Np*K); % get "minus" boundary values

% make edge integration matrix on non-boundary edges
r1D = JacobiGL(0,0,N); V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

nx = reshape(nx,Nfp,Nfaces*K);ny = reshape(ny,Nfp,Nfaces*K);
nxM = nx(:,fpairs(1,:));nyM = ny(:,fpairs(1,:));
nxP = nx(:,fpairs(2,:));nyP = ny(:,fpairs(2,:));

% quiver(x(fM(:)),y(fM(:)),nx(fM(:)),ny(fM(:)))
% plotVerts;quiver(x(vmapM),y(vmapM),nx(:),ny(:))

Efm = sparse(1:NfpT,fM(:),1,NfpT,Np*K);
Efp = sparse(1:NfpT,fP(:),fM~=fP,NfpT,Np*K); % don't subtract off vmapP for boundary nodes, where fM==fP
EfjumpX = spdiag(nxM(:))*Efm + spdiag(nxP(:))*Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms
EfjumpY = spdiag(nyM(:))*Efm + spdiag(nyP(:))*Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms
Ef = [EfjumpX EfjumpY]; % get normal minus jump = sigma_n on face

% get flux ids on boundaries
Ntrp = Nt+1;
% degree of hybrid unknown - Nt = restricted N
if (Nt==0)
    R1D = ones(Nfp,1)/sqrt(Nfp); % interpolate constant to multiple nodes
else    
    V1Dfr = Vandermonde1D(Nt,rtr);
    V1Df = Vandermonde1D(Nt,r1D);                   
    R1D = V1Df/V1Dfr; % interpolate lower order polynomial to order Nf polynomial    
end
sJReduc = reshape(sJ,Nfp,Nfaces*K);
sJReduc = sJReduc(:,fpairs(1,:)); % get unique faces
SJ = spdiag(sJReduc(:));
B = kron(speye(NfacesU),R1D'*M1D)*SJ*Ef; % for fluxes

% return nx,ny to original config
nx = nx(:);
ny = ny(:);
return
