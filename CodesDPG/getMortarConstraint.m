% gets the constraint matrix B for a mortar discretization corresponding to
% the block system
% [Ak B'
%  B  0]
% where Ak is a block diagonal matrix and B enforces orthogonality of the
% solution across faces to polynomials of order Nfr.
% other return values: vmapBF, xf, yf, nxf, nyf = boundary info.
function [B vmapBF xf yf nxf nyf] = getMortarConstraint(Nfr)

Nfr = max(0,Nfr); 

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
Efm = sparse(1:NfpT,fM(:),1,NfpT,Np*K);
Efp = sparse(1:NfpT,fP(:),fM~=fP,NfpT,Np*K); % don't subtract off vmapP for boundary nodes, where fM==fP
Ef = Efm-Efp; % jump matrix (maps volume nodes to jumps over faces) modified to not zero out boundary terms

% make edge integration matrix on non-boundary edges
r1D = JacobiGL(0,0,N); V1D = Vandermonde1D(N,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

% degree of hybrid unknown - Nfr = restricted N
if (Nfr==0)
    R1D = ones(Nfp,1)/sqrt(Nfp); % interpolate constant to multiple nodes
else    
    rfr = JacobiGL(0,0,Nfr);V1Dfr = Vandermonde1D(Nfr,rfr);
    V1Df = Vandermonde1D(Nfr,r1D);                   
    R1D = V1Df*inv(V1Dfr); % interpolate lower order polynomial to order Nf polynomial    
end

sJReduc = reshape(sJ,Nfp,Nfaces*K); 
sJReduc = sJReduc(:,fpairs(1,:)); % get unique faces
SJ = spdiag(sJReduc(:));
B = kron(speye(NfacesU),R1D'*M1D)*SJ*Ef; % for fluxes

% get flux ids on boundaries
Nfrp = Nfr+1;
vmapBF = zeros(Nfrp,NfacesU); vmapBF(:)= 1:Nfrp*NfacesU;
bfaces = any(fM==fP); 
vmapBF = vmapBF(:,bfaces); vmapBF = vmapBF(:);

NfacesB = nnz(bfaces);% get boundary nodes: index into unique face nodes
vmapB = reshape(vmapB,Nfp,NfacesB);
xB = x(vmapB); yB = y(vmapB);
xf = zeros(Nfrp,NfacesB); yf = zeros(Nfrp,NfacesB);
if (Nfr==0)
    xf(:) = (xB(1,:)+xB(Nfp,:))/2;
    yf(:) = (yB(1,:)+yB(Nfp,:))/2;
else % WARNING: only works for straight edge elements. 
    xi = (rfr+1)/2; % interp along 1D points
    for i = 1:Nfrp
        xf(i,:) = xB(1,:) + xi(i)*(xB(Nfp,:)-xB(1,:));
        yf(i,:) = yB(1,:) + xi(i)*(yB(Nfp,:)-yB(1,:));
    end
end
nxf = reshape(nx,Nfp,Nfaces*K); nyf = reshape(ny,Nfp,Nfaces*K);
nxf = nxf(:,fpairs(1,:)); nyf = nyf(:,fpairs(1,:));
nxf = nxf(:,bfaces); nyf = nyf(:,bfaces);

% WARNING: straight edge elem hack - take 1st Nfrp points, all nx/ny same.
nxf = nxf(1:Nfrp,:); nyf = nyf(1:Nfrp,:);

vmapB = vmapB(:);
return
