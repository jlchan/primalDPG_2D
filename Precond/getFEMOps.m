% get P1 matrix (sparsified dense matrix) for nodes on 1 quad element 

function [M1 Dr1 Ds1] = getFEMOps(N)

Globals2D
globals = backupGlobals(1);
[x1 y1] = Nodes2D(N);
tri = delaunay(x1,y1);

% make new mesh from delaunay
EToV = tri;
K = size(tri,1);  Nv = length(x1);
VX = x1(:)';VY = y1(:)';

useQuads = false; Init;  % SWITCH TO TRIANGLES
N = 1; StartUp2D;
[M1, Dr1, Ds1] = getBlockOps();
R = getCGRestriction();

% [r c] = find(R); [ru i] = unique(r);
% xu = x(c(i)); yu = y(c(i));   % get unique nodes

% assemble mass/stiffness matrices
M1 = R*M1*R';
Dr1 = R*Dr1*R';
Ds1 = R*Ds1*R';

useQuads = true; Init % SWITCH BACK TO QUADS
backupGlobals(0, globals);
