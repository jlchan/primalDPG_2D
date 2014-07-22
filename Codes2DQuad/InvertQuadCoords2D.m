function [Rout,Sout] = InvertQuadCoords2D(VX, VY, X, Y, tol)

  % function Rout = InvertQuadCoords2D(VX, VY, X, Y, tol)
% purpose: find coordinates in reference quad of a list of points X in 
%          a physical quad with vertices [X1,X2,X3,X4] 
%          using Newton's method

  Rout = zeros(size(X));
  Sout = zeros(size(X));

  X1 = [VX(1);VY(1)];
  X2 = [VX(2);VY(2)];
  X3 = [VX(3);VY(3)];
  X4 = [VX(4);VY(4)];

  for n=1:length(X)
    Xn = [X(n);Y(n)];

    R = [0;0];

    err = 1;
    
    while(err>tol)

      r = R(1); s = R(2);

      % compute residual for coordinates
      F = 0.25*( (1-r)*(1-s)*X1 + (1+r)*(1-s)*X2 + (1+r)*(1+s)*X3 + (1-r)*(1+s)*X4 ) - Xn;
      
      % Jacobian 
      dFdr = 0.25*( (1-s)*(X2-X1) + (1+s)*(X3-X4) );
      dFds = 0.25*( (1-r)*(X4-X1) + (1+r)*(X3-X2) );
      
      J = [dFdr,dFds];
%      X1,X2,X3,X4,F,Xn,J

      % Newton step
      R = R-(J\F);
      
      % compute error in residual (not updated)
      err = norm(F);
      
    end

    Rout(n) = R(1);
    Sout(n) = R(2);
    
  end

