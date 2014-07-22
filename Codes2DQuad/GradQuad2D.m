function [ux,uy] = Grad2D(u);

% function [ux,uy] = Grad2D(u);
% Purpose: Compute 2D gradient field of scalar u

Globals2D;

% ur = Dr*u; us = Ds*u;
D = full(Dr(1:N+1,1:N+1));

u = reshape(u, N+1,(N+1)*K);

ur = D*u;
ur = reshape(ur, Np, K);

u = transpose(u);
u = reshape(u, N+1, K*(N+1));
us = D*u;

us = reshape(us, K*(N+1), N+1);
us = transpose(us);
us = reshape(us, Np, K);

ux = rx.*ur + sx.*us; uy = ry.*ur + sy.*us;
return
