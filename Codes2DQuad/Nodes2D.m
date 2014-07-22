function [r,s] = Nodes2D(N)

r = JacobiGQ(0, 0, N)*ones(1,N+1);
s = r';
r = r(:); s = s(:);
