% gets the constraint matrix B for a mortar discretization corresponding to
% the block system
% [Ak B'
%  B  0]
% where Ak is a block diagonal matrix and B enforces orthogonality of the
% solution across faces to polynomials of order Nfr.
% other return values: vmapBF, xfB, yfB, nxf, nyf = boundary info.
% xf, yf = face coordinates

function [B vmapBF xf yf nxf nyf fpairs] = getMortarConstraint(Nfr)

Nfr = max(0,Nfr); % cannot have negative flux order

Globals2D

[fM fP fpairs] = getFaceInfo();
NfacesU = size(fpairs,2);

% map field to interface dofs
NfpT = numel(fM); % or fP...
Efm = sparse(1:NfpT,fM(:),1,NfpT,Np*K);
Efp = sparse(1:NfpT,fP(:),fM~=fP,NfpT,Np*K); % don't subtract off vmapP for boundary nodes, where fM==fP
Ef = Efm-Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms

% make edge integration matrix on non-boundary edges
r1D = JacobiGL(0,0,N); V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

% get flux ids on boundaries
[rfr xf yf nxf nyf] = getFaceNodes(Nfr,fM,fpairs);
Nfrp = Nfr+1;

% degree of hybrid unknown - Nfr = restricted N
if (Nfr==0)
    R1D = ones(Nfp,1)/sqrt(Nfp); % interpolate constant to multiple nodes
else    
    V1Dfr = Vandermonde1D(Nfr,rfr);
    V1Df = Vandermonde1D(Nfr,r1D);                   
    R1D = V1Df/V1Dfr; % interpolate lower order polynomial to order Nf polynomial    
end

sJReduc = reshape(sJ,Nfp,Nfaces*K);
sJReduc = sJReduc(:,fpairs(1,:)); % get unique faces
SJ = spdiag(sJReduc(:));
B = kron(speye(NfacesU),R1D'*M1D)*SJ*Ef; % for fluxes

% boundary nodes
vmapBF = zeros(Nfrp,NfacesU); vmapBF(:)= 1:Nfrp*NfacesU;
bfaces = any(fM==fP); 
vmapBF = vmapBF(:,bfaces); vmapBF = vmapBF(:);


return
