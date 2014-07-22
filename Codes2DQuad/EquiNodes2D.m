function [r,s] = EquiNodes2D(N)

r = linspace(-1, 1, N+1);
r = r(:)*ones(1,N+1);
s = r';
r = r(:); s = s(:);
