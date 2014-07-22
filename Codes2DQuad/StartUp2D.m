% Purpose : Setup script, building operators, grid, metric,
%           and connectivity tables.
warning('Quad version being run')

N = N;

% Definition of constants
Nfp = N+1; Np = (N+1)*(N+1); 
Nfaces=4; NODETOL = 1e-12;

% Compute nodal set
r = JacobiGL(0, 0, N)*ones(1,N+1);
s = r';
r = r(:); s = s(:);

% Build reference element matrices
V = Vandermonde2D(N,r,s);
invV = inv(V);

MassMatrix = invV'*invV;
[Dr,Ds] = Dmatrices2D(N, r, s, V);

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; 
vc = EToV(:,3)'; vd = EToV(:,4)'; 

x = 0.25*((1-r).*(1-s)*VX(va)+...
	   (1+r).*(1-s)*VX(vb)+...
	   (1+r).*(1+s)*VX(vc)+...
	   (1-r).*(1+s)*VX(vd));

y = 0.25*((1-r).*(1-s)*VY(va)+...
	   (1+r).*(1-s)*VY(vb)+...
	   (1+r).*(1+s)*VY(vc)+...
	   (1-r).*(1+s)*VY(vd));

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r-1) < NODETOL)';
fmask3   = find( abs(s-1) < NODETOL)';
fmask4   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

Fx = x(Fmask(:), :);
Fy = y(Fmask(:), :);

% Create surface integral terms
LIFT = LiftQuad2D();

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = NormalsQuad2D();
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = tiConnectQuad2D(EToV);

% Build connectivity maps
BuildMapsQuad2D;

return
