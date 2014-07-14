% [rfr xf yf nxf nyf] = getFaceNodes(Nfr,fM,fpairs)
% rfr = reference face nodes in 2D [-1,1]
% xf,yf = face node positions
% nxf,nyf = face node normals (from interior/minus perspective)

function [rfr xf yf nxf nyf] = getFaceNodes(Nfr,fM,fpairs)

Globals2D

% reference face nodes

NfacesU = size(fpairs,2);
Nfrp = Nfr+1;

xf = zeros(Nfrp,NfacesU); yf = zeros(Nfrp,NfacesU);
if (Nfr==0)
    rfr = 1; 
    xf(:) = (x(fM(1,:))+x(fM(Nfp,:)))/2;
    yf(:) = (y(fM(1,:))+y(fM(Nfp,:)))/2;    
else % WARNING: only works for straight edge elements. 
    rfr = JacobiGL(0,0,Nfr);
    xi = (rfr+1)/2; % interp along 1D points
    for i = 1:Nfrp
        xf(i,:) = x(fM(1,:)) + xi(i)*(x(fM(Nfp,:))-x(fM(1,:)));
        yf(i,:) = y(fM(1,:)) + xi(i)*(y(fM(Nfp,:))-y(fM(1,:)));        
    end
end
nxf = reshape(nx,Nfp,Nfaces*K); nyf = reshape(ny,Nfp,Nfaces*K);
nxf = nxf(:,fpairs(1,:)); nyf = nyf(:,fpairs(1,:));
nxf = nxf(1:Nfrp,:); nyf = nyf(1:Nfrp,:); % WARNING: straight edge elem hack - take 1st Nfrp points, all nx/ny same.
