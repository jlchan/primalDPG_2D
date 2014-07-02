function plotNodes()
Globals2D
figure
plot(x,y,'.');hold on;
j = 1;
for k = 1:K
    for i = 1:size(x,1)
        off = .2*((k-1)/K);
        text(x(i,k)+off,y(i,k)+off,num2str(j),'fontsize',16)
        j = j+1;
    end
end
% PlotMesh2D

