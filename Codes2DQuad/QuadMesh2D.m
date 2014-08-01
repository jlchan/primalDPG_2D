% builds a simple mesh of NK-by-NK quad elements 

function [Nv, VX, VY, K, EToV] = QuadMesh2D(NK)

Nx = NK;Ny = NK;

Nxp = Nx+1;
Nyp = Ny+1;
Nv = Nxp*Nyp;
K = Nx*Ny;

x1D = linspace(-1,1,Nxp);
[y, x] = meshgrid(x1D);
[I, J] = meshgrid(1:Nxp);
inds = (I-1)*Ny + (J+I-1);
EToV = zeros(K,4);
k = 1;
for i = 1:Nx
    for j = 1:Ny
        EToV(k,:) = [inds(i,j) inds(i+1,j) inds(i+1,j+1) inds(i,j+1)  ];
        k = k+1;
    end
end

VX = x(:)';
VY = y(:)';


