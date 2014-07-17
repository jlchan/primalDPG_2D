Globals2D;
FaceGlobals2D;

% N = 5; Nf = 5;  % when N = even, Nf = N-1, fails?
% % [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% %     Nv = 3;  VX = VX(EToV(1,:)); VY = VY(EToV(1,:));
% %     EToV = [3 1 2];  K = 1;
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% StartUp2D;

if (Nf>N) 
    warning('Nf > N.')
end

Nfrp = Nf+1;

[fM fP fpairs] = getFaceInfo();
NfacesU = size(fpairs,2);

% get flux ids on boundaries
[rfr xf yf nxf nyf] = getFaceNodes(Nf,fM,fpairs);

% fmap = Nfaces*Nfrp x K array of face node dofs
sharedFaces = ~ismember(fpairs(2,:),fpairs(1,:));
fmap = zeros(Nfrp,Nfaces*K);
fmap(:,fpairs(1,:)) = reshape(1:Nfrp*NfacesU,Nfrp,NfacesU);
fmap(:,fpairs(2,sharedFaces)) = fmap(:,fpairs(1,sharedFaces));
fmap = reshape(fmap,Nfrp*Nfaces,K);

fmapB = reshape(1:Nfrp*NfacesU,Nfrp,NfacesU); 
bfaces = any(fM==fP); 
fmapB = fmapB(:,bfaces); 
fmapB = fmapB(:);

% plotVerts
% % plot(x,y,'rs')
% % plot(x(vmapB),y(vmapB),'k*','markersize',12)
% plot(xf,yf,'g^')
% plot(xf(fmapB),yf(fmapB),'mo','markersize',12)