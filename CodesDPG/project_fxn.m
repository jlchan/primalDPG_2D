function f = project_fxn(fhandle,Corder)

Globals2D;

% compute projections of exact solutions for bdata
[cubR,cubS,cubW, Ncub] = Cubature2D(Corder); Vcub = Vandermonde2D(N,cubR,cubS);
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
xcub = 0.5*(-(cubR+cubS)*VX(va)+(1+cubR)*VX(vb)+(1+cubS)*VX(vc));
ycub = 0.5*(-(cubR+cubS)*VY(va)+(1+cubR)*VY(vb)+(1+cubS)*VY(vc));

Interp = Vcub*invV; % interp to cubature points
Wcub = diag(cubW);
Minv = V*V';
% uex_f = @(x,y) sin(pi*x).*sin(pi*y);
f = Minv*Interp'*Wcub*fhandle(xcub,ycub);
f = f(:);
