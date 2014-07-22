function [Nv, VX, VY, K, EToV] = MeshToQuad2D()

Globals2D;
Globals2D;

CX = sum(VX(EToV), 2)/Nfaces;
CY = sum(VY(EToV), 2)/Nfaces;
Nnodes = length(VX);

VX = [VX,CX'];
VY = [VY,CY'];

vnum = [ 1, 2; 2, 3; 3, 1];

K = 0;
EToV = zeros(K, 4);
for k1=1:K
  for f1=1:Nfaces
    k2 = EToE(k1,f1);
    f2 = EToF(k1,f1);
    if(k1<k2)
      K = K+1;

      v1 = Nnodes+k1;
      v2 = EToV(k1,vnum(f1,1));
      v3 = Nnodes+k2;
      v4 = EToV(k1,vnum(f1,2));
      
      EToV(K,:) = [v1,v2,v3,v4];
    end
  end
end

Nv = length(VX);
