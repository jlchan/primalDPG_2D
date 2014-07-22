function refc = BuildIntegrationQuad2D(intC, intG)

Globals2D;
Globals2D;

% 1.1 Extract cubature nodes and weights
[gz,gw] = JacobiGQ(0, 0, intC);
Ng = length(gz);
cR = gz*ones(1,Ng);
cS = cR';
cW = gw*gw';
Ncub = length(cW(:));
cR = cR(:); 
cS = cS(:);
cW = cW(:);

% 1.1. Build interpolation matrix (nodes->cubature nodes)
cV = InterpMatrixuad2D(cR, cS);
  
% 1.2 Evaluate derivatives of Lagrange interpolants at cubature nodes
[cDr,cDs] = Dmatrices2D(N, cR, cS, V);

% differentiate coordinates
xr = Dr*x; yr = Dr*y;
xs = Ds*x; ys = Ds*y;

cI = InterpMatrixuad2D(cR, cS);

% evaluate these derivatives at cubature 
cxr = cI*xr; cxs = cI*xs;
cyr = cI*yr; cys = cI*ys;

% differentiate again with results evaluated at cubature
cxrr = cDr*xr; cxrs = cDr*xs; cxsr = cDs*xr; cxss = cDs*xs;
cyrr = cDr*yr; cyrs = cDr*ys; cysr = cDs*yr; cyss = cDs*ys;

% construct geometric factors
cJ = cxr.*cys - cxs.*cyr;
crx =  cys./cJ; csx = -cyr./cJ;
cry = -cxs./cJ; csy =  cxr./cJ;

csqrtJ = sqrt(cJ);

% now manually differentiate Jacobian
cdJdr = cxrr.*cys + cxr.*cysr - cxsr.*cyr - cxs.*cyrr;
cdJds = cxrs.*cys + cxr.*cyss - cxss.*cyr - cxs.*cyrs;

% chain rule
cdJdx = crx.*cdJdr + csx.*cdJds;
cdJdy = cry.*cdJdr + csy.*cdJds;

cdlogJdx = cdJdx./cJ; 
cdlogJdy = cdJdy./cJ; 

cV = InterpMatrixuad2D(cR, cS);

refc.V = cV; % interpolates nodes to cubature
refc.Dr = cDr; % differentiate nodal data and evaluate at cubature nodes
refc.Ds = cDs; % differentiate nodal data and evaluate at cubature nodes

refc.DrT = V*V'*cDr'*diag(cW);
refc.DsT = V*V'*cDs'*diag(cW);

refc.Proj = V*V'*cV'*diag(cW); % project from cubature data to nodes (no Jacobians built in)
refc.dlogJdx = cdlogJdx;
refc.dlogJdy = cdlogJdy;
refc.drdx = crx;
refc.drdy = cry;
refc.dsdx = csx;
refc.dsdy = csy;
refc.J    = cJ;

refc.W = cW;
