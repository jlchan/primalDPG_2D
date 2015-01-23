Globals2D;
N = 4;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');
StartUp2D;

%%
PlotMesh2D()
hold on
% no overlap/face overlap only
for k = K-2
    elems = unique([k EToE(k,:)]); % face overlap only
    %     inds = unique(rr(:,elems));
    for e = 1:length(elems)
        vi = [EToV(elems(e),:) EToV(elems(e),1)];
        plot(VX(vi),VY(vi),'r-','linewidth',3)
    end    
end
print(gcf,'-depsc','../figs/facePatch.eps')

%%
clf
PlotMesh2D()
hold on

[R, vmapBT] = getCGRestriction();
Norder = N;
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
for i = length(vnodes)-1    
    [~, elems] = find(vnodes(i)==vertices);
    inds = unique(rr(:,elems));
    for e = 1:length(elems)
        vi = [EToV(elems(e),:) EToV(elems(e),1)];
        plot(VX(vi),VY(vi),'r-','linewidth',3)
    end
end
print(gcf,'-depsc','../figs/vertexPatch.eps')
