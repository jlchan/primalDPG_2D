
function [test,xint] = EdgeIntersect3D(a1, b1, a2, b2)

  %
  %  a1 + (b1-a1)*t1 = a2 + (b2-a2)*t2
  %
  %  (b2-a2)x(b1-a1)*t1 = (b2-a2)x(a2-a1)
  %
  %  t1 = ((b2-a2)x(b1-a1)).((b2-a2)x(a2-a1)) / (||(b2-a2)x(b1-a1)||^2)
  %  = c.(-d)/(c.c)
  %
  %  (b1-a1)x(b2-a2)*t2 = (b1-a1)x(a1-a2)
  %
  %  t2 = ((b1-a1)x(b2-a2)).((b1-a1)x(a1-a2)) / (||(b2-a2)x(b1-a1)||^2)
  %  = (-c).(-e)/(c.c)
  %

  test = false;
  xint = -999;

  cx  = (b2(2)-a2(2))*(b1(3)-a1(3))-(b2(3)-a2(3))*(b1(2)-a1(2));
  cy  = (b2(3)-a2(3))*(b1(1)-a1(1))-(b2(1)-a2(1))*(b1(3)-a1(3));
  cz  = (b2(1)-a2(1))*(b1(2)-a1(2))-(b2(2)-a2(2))*(b1(1)-a1(1));
  
  dx  = (b2(2)-a2(2))*(a1(3)-a2(3))-(b2(3)-a2(3))*(a1(2)-a2(2));
  dy  = (b2(3)-a2(3))*(a1(1)-a2(1))-(b2(1)-a2(1))*(a1(3)-a2(3));
  dz  = (b2(1)-a2(1))*(a1(2)-a2(2))-(b2(2)-a2(2))*(a1(1)-a2(1));
  
  ex  = (b1(2)-a1(2))*(a2(3)-a1(3))-(b1(3)-a1(3))*(a2(2)-a1(2));
  ey  = (b1(3)-a1(3))*(a2(1)-a1(1))-(b1(1)-a1(1))*(a2(3)-a1(3));
  ez  = (b1(1)-a1(1))*(a2(2)-a1(2))-(b1(2)-a1(2))*(a2(1)-a1(1));
  
  t1  = -(cx*dx+cy*dy+cz*dz);
  t2  = cx*ex+cy*ey+cz*ez;
  
  mag = cx*cx+cy*cy+cz*cz;
  
  umTOL4 = 1e-6;
  if(mag < umTOL4)
    % lines are colinear --
    %   should be caught by "in test" 

    return;
  end

  t1 = t1/mag;
  t2 = t2/mag;

  x1 = zeros(3,1); 
  x2 = zeros(3,1);

  x1(1) = a1(1)+t1*(b1(1)-a1(1));
  x1(2) = a1(2)+t1*(b1(2)-a1(2));
  x1(3) = a1(3)+t1*(b1(3)-a1(3));

  x2(1) = a2(1)+t2*(b2(1)-a2(1));
  x2(2) = a2(2)+t2*(b2(2)-a2(2));
  x2(3) = a2(3)+t2*(b2(3)-a2(3));

  % make sure the segments really intersect.
  d = (x1(1)-x2(1))^2+(x1(2)-x2(2))^2+(x1(3)-x2(3))^2;

  umTOL4 = 1e-6;

  if(d>umTOL4)
    return;
  end
  
  if( (t1 > -umTOL4)  & (t2 > -umTOL4) & (t1 < 1.0+umTOL4) & (t2 < 1.0+umTOL4))
    % woo-hoo we have intersection

    xint = x1;
    test = true;
    return;
  end
  
  return;
