% get P1 matrix (sparsified dense matrix) for nodes 

function [M1 Dx1 Dy1] = getFEMOps(x1,y1)

Globals2D
globals = backupGlobals(1);

tri = delaunay(x1,y1);

% make new mesh from delaunay
EToV = tri;
K = size(tri,1);  Nv = length(x1);
VX = x1(:)';VY = y1(:)';

N = 1; StartUp2D;
[M, Dx, Dy] = getBlockOps();
R = getCGRestriction();
% [r c] = find(R); [ru i] = unique(r);
% xu = x(c(i)); yu = y(c(i));   % get unique nodes

% assemble mass/stiffness matrices
M1 = R*M*R';
Dx1 = R*Dx*R';
Dy1 = R*Dy*R';

backupGlobals(0, globals);
