function [Nv, VX, VY, K, EToV] = TriToQuad2D(FileName)

Nfaces = 3;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(FileName);
[EToE,EToF]= tiConnect2D(EToV);

Nfaces = 4;

% split all elements
newEToV = zeros(3*K,Nfaces);
newK = 1;
newVX = VX;
newVY = VY;
for k=1:K
  
  va = EToV(k,1);  vb = EToV(k,2);  vc = EToV(k,3);
  ka = EToE(k,1);  kb = EToE(k,2);  kc = EToE(k,3);
  fa = EToF(k,1);  fb = EToF(k,2);  fc = EToF(k,3);

  vab = Nv + 1 + max(k*(Nfaces+1) + 1, ka*(Nfaces+1) + fa);
  vbc = Nv + 1 + max(k*(Nfaces+1) + 2, kb*(Nfaces+1) + fb);
  vca = Nv + 1 + max(k*(Nfaces+1) + 3, kc*(Nfaces+1) + fc);

  vcent = Nv + 1 + (Nfaces+1)*K + 1  + k; 

  newEToV(newK, :) = [va,vab,vcent,vca]; newK = newK + 1;
  newEToV(newK, :) = [vab,vb,vbc,vcent]; newK = newK + 1;
  newEToV(newK, :) = [vca,vcent,vbc,vc]; newK = newK + 1;
  
  xa = VX(va);  xb = VX(vb);   xc = VX(vc); 
  ya = VY(va);  yb = VY(vb);   yc = VY(vc); 

  newVX(vab) = (xa+xb)/2;
  newVX(vbc) = (xb+xc)/2;
  newVX(vca) = (xc+xa)/2;

  newVY(vab) = (ya+yb)/2;
  newVY(vbc) = (yb+yc)/2;
  newVY(vca) = (yc+ya)/2;

  newVX(vcent) = (xa+xb+xc)/3;
  newVY(vcent) = (ya+yb+yc)/3;
	
end

ids = sort(unique(newEToV(:)), 'ascend');

map = zeros(max(ids));
map(ids) = 1:length(ids);

EToV = map(newEToV);
VX(map(ids)) = newVX(ids);
VY(map(ids)) = newVY(ids);
K = newK-1;
Nv = length(VX);
