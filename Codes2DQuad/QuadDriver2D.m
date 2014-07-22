% Driver script for solving the 2D vacuum Maxwell's equations on TM form
Globals2D;
Globals2D;

% Polynomial order used for approximation 
N = 6;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% build Quad mesh
[Nv, VX, VY, K, EToV] = MeshToQuad2D();

StartUpQuad2D;

EQ = sin(pi*x).*sin(pi*y);
PlotFieldQuad2D(N, x, y, EQ);
