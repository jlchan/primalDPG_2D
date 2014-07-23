function [r,s] = Nodes2D(N)

r = JacobiGL(0, 0, N)*ones(1,N+1);
s = r';
r = r(:); s = s(:);
