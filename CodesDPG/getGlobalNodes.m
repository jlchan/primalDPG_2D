function [xg yg] = getGlobalNodes(rl,sl)

Globals2D

quads = size(EToV,2)==4;
if quads
    va = EToV(:,1)'; vb = EToV(:,2)';
    vc = EToV(:,3)'; vd = EToV(:,4)';
    
    xg = 0.25*((1-rl).*(1-sl)*VX(va)+...
        (1+rl).*(1-sl)*VX(vb)+...
        (1+rl).*(1+sl)*VX(vc)+...
        (1-rl).*(1+sl)*VX(vd));
    
    yg = 0.25*((1-rl).*(1-sl)*VY(va)+...
        (1+rl).*(1-sl)*VY(vb)+...
        (1+rl).*(1+sl)*VY(vc)+...
        (1-rl).*(1+sl)*VY(vd));
else % triangles
    va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
    xg = 0.5*(-(rl+sl)*VX(va)+(1+rl)*VX(vb)+(1+sl)*VX(vc));
    yg = 0.5*(-(rl+sl)*VY(va)+(1+rl)*VY(vb)+(1+sl)*VY(vc));
end