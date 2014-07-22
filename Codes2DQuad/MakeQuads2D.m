function [Nv, VX, VY, K, EToV] = MakeQuads2D(Nx)

x1d = linspace(-1, 1, Nx+1);

VX = zeros(Nx+1,Nx+1);
VY = zeros(Nx+1,Nx+1);
EToV = zeros( Nx*Nx, 4);
sk = 1;

for n=1:Nx+1
  for m=1:Nx+1
    VX(n,m) = x1d(n);
    VY(n,m) = x1d(m);

  end
end

for n=2:2:Nx
	sg = -1; fac = 1/Nx;
  for m=1:1:Nx+1
    VX(n,m) = VX(n,m) + sg*fac;
    sg = -sg;
  end
end

%VY(:,2) = VY(:,2) + .25*[-1;1;-1;1;-1];
%VY(:,4) = VY(:,4) + .25*[-1;1;-1;1;-1];

for n=1:Nx
  for m=1:Nx
    EToV(sk,1) = (n  )+(m-1)*(Nx+1);
    EToV(sk,2) = (n+1)+(m-1)*(Nx+1);
    EToV(sk,3) = (n+1)+(m)*(Nx+1);
    EToV(sk,4) = (n  )+(m)*(Nx+1); 
    
    sk = sk+1;

  end
end

VX = VX(:)';
VY = VY(:)';

K = sk-1;
Nv = (Nx+1)*(Nx+1);
