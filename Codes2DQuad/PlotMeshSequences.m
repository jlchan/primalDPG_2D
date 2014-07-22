% Driver script for solving the 2D vacuum Maxwell's equations on TM form
% Globals2D;
Globals2D;
N = 1;


for good=1
    
[Nv, VX, VY, K, EToV] = MakeQuads2D(2);

for loop=1:4 
StartUpQuad2D;

if(good)
if(loop>1)
CheapRefine2D;
end
else
[Nv, VX, VY, K, EToV] = MakeQuads2D(2^(loop+1));
end

StartUpQuad2D;

PlotMeshQuad2D();
axis equal; axis off; set(gcf, 'Color', 'White')
drawnow; pause(1)
hold on
ids = EToV(1,:);
fill(VX(ids),VY(ids), 'b')
hold off

% print('-dpdf', '-painters', sprintf('mesh_%d_%d.pdf', good, loop));

end
end
