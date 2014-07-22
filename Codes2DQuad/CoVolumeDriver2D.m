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

% build transfer matrices
if(1==1)
[TToQ, QToT] = BuildProjection2D();
else
TToQ = [];
QToT = [];
end

% quick test of Maxwell (yes - dubious basis )
FinalTime = 1;

H = 0*rand(Np, K, 2);
E = exp(-10*(x.^2+y.^2)); % +.1*rand(Np, K);
%E = sin(pi*x).*sin(pi*y);
%E = ones(Np, K);
% [H, E] = Maxwell2D(H, E, FinalTime, TToQ, QToT);

F = BuildOverlap2D();
% [H,E] = MaxwellOverlap2D(H, E, FinalTime, F);

HQ = 0*rand(Np, K, 2);
EQ = exp(-30*(x.^2+y.^2)); % +.1*rand(Np, K);
[HQ, EQ] = MaxwellQuad2D(HQ, EQ, FinalTime, TToQ, QToT);

%max(max(abs(E(:)-QToT*EQ(:))))

cpt = cputime;
%[HQ, E] = MaxwellMixed2D(HQ, E, FinalTime, TToQ, QToT);
mixedtime = cputime-cpt
