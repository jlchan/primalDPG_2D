% get P1 matrix (sparsified dense matrix) for nodes 


function [M1 Dx1 Dy1] = getFEMOps(N)

Globals2D
globals = backupGlobals(1);

[x1 y1] = Nodes2D(N);
tri = delaunay(x1,y1);
% REMOVE SIZE ZERO TRIANGLES - this code doesn't work yet

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

% % scale with geom factors
% blkDr = kron(speye(K),Dr);
% blkDs = kron(speye(K),Ds);
% blkM = kron(speye(K),MassMatrix);
% 
% M = spdiag(J(:))*blkM; % J = h^2
% Dx = spdiag(rx(:))*blkDr + spdiag(sx(:))*blkDs;
% Dy = spdiag(ry(:))*blkDr + spdiag(sy(:))*blkDs;
