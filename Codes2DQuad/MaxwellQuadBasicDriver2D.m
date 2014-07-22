% Driver script for solving the 2D vacuum Maxwell's equations on TM form
% Globals2D;
Globals2D;

% Polynomial order used for approximation 
N = 4;

% % Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
% 
% % Initialize solver and construct grid and metric
% StartUp2D;
% 
% % build Quad mesh
% [Nv, VX, VY, K, EToV] = MeshToQuad2D();

% [Nv, VX, VY, K, EToV] = MakeQuads2D(1);
[Nv, VX, VY, K, EToV] = QuadMesh2D(16);
StartUp2D;

EQ = sin(pi*x).*sin(pi*y);
PlotFieldQuad2D(N, x, y, EQ);

% quick test of Maxwell (yes - dubious basis )
FinalTime = 1;

HQ = 0*rand(Np, K, 2);
EQ = exp(-30*(x.^2+y.^2));% +.1*rand(Np, K);
[HQ, EQ] = MaxwellQuadBasic2D(HQ, EQ, FinalTime);

