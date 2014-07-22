function [r,s] = InvertTriCoords(VX, VY, X, Y)

% function [r,s] = InvertTriCoords(VX, VY, X, Y)
% purpose: find local (r,s) coordinates in the k'th element of given coordinates
%          [only works for straight sided triangles]

X1 = [VX(1);VY(1)];
X2 = [VX(2);VY(2)];
X3 = [VX(3);VY(3)];

A = [ X2-X1, X3-X1];

M = length(X);

R = A\(2*[X';Y'] - (X2+X3)*ones(1,M) );

r = R(1,:)';
s = R(2,:)';

return
