% plane waves

x = linspace(-1,1,150);
[X Y] = meshgrid(x);
k = 25;
t = pi/4;
U = exp(i*k*(cos(t)*X+sin(t)*Y));

% t = pi/3;
% U = U + exp(i*k*(cos(t)*X+sin(t)*Y));
% t = pi/2;
% U = U + exp(i*k*(cos(t)*X+sin(t)*Y));
% t = pi;
% U = U + exp(i*k*(cos(t)*X+sin(t)*Y));

imagesc(real(U))
shading interp