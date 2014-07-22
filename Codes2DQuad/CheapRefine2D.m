
% split all elements
Nv = length(VX);
newEToV = zeros(4*K,Nfaces);
newK = 1;
newVX = VX;
newVY = VY;
for k=1:K
  
  va = EToV(k,1);
  vb = EToV(k,2);
  vc = EToV(k,3);
  vd = EToV(k,4);

  ka = EToE(k,1);
  kb = EToE(k,2);
  kc = EToE(k,3);
  kd = EToE(k,4);

  fa = EToF(k,1);
  fb = EToF(k,2);
  fc = EToF(k,3);
  fd = EToF(k,4);

  vab = Nv + 1 + max(k*(Nfaces+1) + 1, ka*(Nfaces+1) + fa);
  vbc = Nv + 1 + max(k*(Nfaces+1) + 2, kb*(Nfaces+1) + fb);
  vcd = Nv + 1 + max(k*(Nfaces+1) + 3, kc*(Nfaces+1) + fc);
  vda = Nv + 1 + max(k*(Nfaces+1) + 4, kd*(Nfaces+1) + fd);

  vcent = Nv + 1 + (Nfaces+2)*K + 1  + k; 

  newEToV(newK, :) = [va,vab,vcent,vda]; newK = newK + 1;
  newEToV(newK, :) = [vab,vb,vbc,vcent]; newK = newK + 1;
  newEToV(newK, :) = [vcent,vbc,vc,vcd]; newK = newK + 1;
  newEToV(newK, :) = [vda,vcent,vcd,vd]; newK = newK + 1;
  
  xa = VX(va);  xb = VX(vb);   xc = VX(vc);   xd = VX(vd);
  ya = VY(va);  yb = VY(vb);   yc = VY(vc);   yd = VY(vd);

  newVX(vab) = (xa+xb)/2;
  newVX(vbc) = (xb+xc)/2;
  newVX(vcd) = (xc+xd)/2;
  newVX(vda) = (xd+xa)/2;

  newVY(vab) = (ya+yb)/2;
  newVY(vbc) = (yb+yc)/2;
  newVY(vcd) = (yc+yd)/2;
  newVY(vda) = (yd+ya)/2;

  newVX(vcent) = (xa+xb+xc+xd)/4;
  newVY(vcent) = (ya+yb+yc+yd)/4;
	
end

ids = sort(unique(newEToV(:)), 'ascend');

map = zeros(max(ids));
map(ids) = 1:length(ids);

EToV = map(newEToV);
VX(map(ids)) = newVX(ids);
VY(map(ids)) = newVY(ids);
K = newK-1;
