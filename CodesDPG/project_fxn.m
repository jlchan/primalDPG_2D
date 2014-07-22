function f = project_fxn(fhandle,Corder)

Globals2D;

% compute projections of exact solutions for bdata
[cubR,cubS,cubW, Ncub] = Cubature2D(Corder); Vcub = Vandermonde2D(N,cubR,cubS);

[xcub ycub] = getGlobalNodes(cubR,cubS);

Interp = Vcub*invV; % interp to cubature points
Wcub = diag(cubW);
Minv = V*V';
% uex_f = @(x,y) sin(pi*x).*sin(pi*y);
f = Minv*Interp'*Wcub*fhandle(xcub,ycub);
f = f(:);
