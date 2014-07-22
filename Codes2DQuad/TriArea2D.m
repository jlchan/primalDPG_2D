
function A = TriArea2D(x, y)

  % function A = TriArea2D(x, y)
  % purpose: compute area of triangle given in triples x,y
     
  x12 = x(2)-x(1);
  y12 = y(2)-y(1);
  x13 = x(3)-x(1);
  y13 = y(3)-y(1);

  A = 0.5*(x12*y13-y12*x13);

  return;
