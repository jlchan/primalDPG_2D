
function [outflag, alpha, beta, gam] = FindIntersectionLine2D(vx, vy, sx, sy)

  outflag = 0;
  zertol = 0;
  onetol = 1;
  Nfaces = 3;
  
  x1 = [sx(1);sy(1)];
  x2 = [sx(2);sy(2)];
  
  flag = zeros(Nfaces,1);
  flagXYZ = zeros(2, Nfaces);
  vnum = [1,2;2,3;3,1];
  for f=1:Nfaces
    
    id1 = vnum(f,1);
    id2 = vnum(f,2);
    
    v1 = [ vx(id1); vy(id1)];
    v2 = [ vx(id2); vy(id2)];
    
    % intersection test: http://www.Graphics.Cornell.EDU/pubs/1997/MT97.pdf
    % point T(u,v) on a triangle is given by
    % T(u,v) = (1-u-v)*v1 + u*v2 + v*v3
    % intersection:   x1 + t*(x2-x1) = (1-u-v)*v1 + u*v2 + v*v3
    % [ -(x2-x1), v2-v1, v3-v1] [ t; u; v] = x1 - v1
    
    tuv = [-(x2-x1), v2-v1]\(x1-v1);
    
    if(min(tuv)>=zertol && max(tuv)<=onetol)
      % yeah -- intersection add this segment to list on k1
      flag(f) = 1;
      flagXYZ(:,f) = (1-tuv(2))*v1 + tuv(2)*v2;
    end
  end   

  
  alpha = 0; beta = 0; gam = 0; 
  if(sum(flag)==2)
    ids = find(flag);

    ix1 = flagXYZ(1,ids(1));
    iy1 = flagXYZ(2,ids(1));
    ix2 = flagXYZ(1,ids(2));
    iy2 = flagXYZ(2,ids(2));
    alpha = -(iy2-iy1);
    beta  = ix2-ix1;
    gam = iy2*ix1-ix2*iy1;

alpha*ix1+beta*iy1+gam;
alpha*ix2+beta*iy2+gam;

    outflag = 1;
  else
    outflag = 0;
  end

  if(abs(alpha)+abs(beta)<1e-10)
    outflag = 0;
  end
