function plotVerts()

Globals2D
% figure
plot(VX,VY,'.');hold on;
Nv = length(VX);
for i = 1:Nv
    off = .1*((i-1)/Nv);
    text(VX(i)+off,VY(i)+off,num2str(i),'fontsize',16,'color','r')        
end

K = size(EToV,1);
for k = 1:K
    cx = sum(VX(EToV(k,:)))/3;
    cy = sum(VY(EToV(k,:)))/3;
    text(cx-1/K,cy,['K=' num2str(k)],'fontsize',16)
    v = EToV(k,:);
    for j = 1:2        
        plot(VX([v(j) v(j+1)]),VY([v(j) v(j+1)]))
    end
    plot(VX([v(3) v(1)]),VY([v(3) v(1)]))
end
% PlotMesh2D

