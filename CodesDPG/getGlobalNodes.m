function [xg yg] = getGlobalNodes(rl,sl)

Globals2D

quads = 0;
if quads    
    error('Nothing implemented for quads yet')
else % triangles    
    va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
    xg = 0.5*(-(rl+sl)*VX(va)+(1+rl)*VX(vb)+(1+sl)*VX(vc));
    yg = 0.5*(-(rl+sl)*VY(va)+(1+rl)*VY(vb)+(1+sl)*VY(vc));
end