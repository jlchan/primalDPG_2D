function [rhsHQ, rhsEQ] = MaxwellQuadBasicRHS2D(refc, HQ,EQ)

  % function [rhsHQ, rhsEQ] = MaxwellRHS2D(refc, HQ,EQ)
% Purpose  : Evaluate RHS flux in 2D Maxwell TM form 

Globals2D;

% approximate with H <- sqrt(1/J)*Hm, E <- sqrt(1/J)*Em, test <- sqrt(1/J)*test
% solve for Hm ~ H*sqrt(J)

%i.   (phi/sqrt(J), d/dt (H/sqrt(J)) J)_Dhat = 
%ii.  (phi/sqrt(J), -curl (E/sqrt(J)) J)_Dhat 
%iii. (phi/sqrt(J), -n x [(E/sqrt(J))] sJ)_Dhat ...

cV = refc.V; cDr = refc.Dr; cDs = refc.Ds; cProj = refc.Proj;

cHx = cV*Hx;
cHy = cV*Hy;
cEz = cV*Ez;

cdHxdr = cDr*Hx; cdHxds = cDs*Hx;
cdHydr = cDr*Hy; cdHyds = cDs*Hy;
cdEzdr = cDr*Ez; cdEzds = cDs*Ez; 

% i. (phi, d/dt (H) )_Dhat = 

% ii.  (phi, -curl E)_Dhat + (phi/sqrt(J), -grad (1/sqrt(J)) x E J)_Dhat
%    = (phi, -curl E)_Dhat + (phi/sqrt(J), 0.5*grad J/J^(3/2) x E J)_Dhat
%    = (phi, -curl E)_Dhat + (phi, 0.5*grad J/J x E )_Dhat

% no Jacobian (cubature to handle rational geo factors)
rhsHx = cProj*( -0.5*refc.dlogJdy.*cEz ) + refc.DrT*(refc.drdy.*cEz) + refc.DsT*(refc.dsdy.*cEz);
rhsHy = cProj*(  0.5*refc.dlogJdx.*cEz ) - refc.DrT*(refc.drdx.*cEz) - refc.DsT*(refc.dsdx.*cEz);

rhsEz =  cProj*(  0.5*refc.dlogJdy.*cHx  - refc.drdy.*cdHxdr - refc.dsdy.*cdHxds);
rhsEz =  cProj*( -0.5*refc.dlogJdx.*cHy  + refc.drdx.*cdHydr + refc.dsdx.*cdHyds) + rhsEz;

% 
gaussI = gauss.interp;
gJ = gauss.J;
gsJ = gauss.sJ;

gmapP = gauss.mapP;
gmapM = gauss.mapM;
gmapB = gauss.mapB;

gnx = gauss.nx;
gny = gauss.ny;

gHx = gaussI*Hx;
gHy = gaussI*Hy;
gEz = gaussI*Ez;

% iii. (phi/sqrt(J), -n x [(E/sqrt(J))] sJ)_Dhat ...
gsqJ = sqrt(gJ);
gHx = gHx./gsqJ;
gHy = gHy./gsqJ;
gEz = gEz./gsqJ;

gdHx = gHx(gmapM)-gHx(gmapP);
gdHy = gHy(gmapM)-gHy(gmapP);
gdEz = gEz(gmapM)-gEz(gmapP);

gaEz = -gEz(gmapM)-gEz(gmapP);

% make change here for exact solution (check isempty HxBC)
gdHx(gmapB) = 0;
gdHy(gmapB) = 0;
gdEz(gmapB) = 2*gEz(gmapB);
gaEz(gmapB) = 0;

gndotdH =  gnx.*gdHx+gny.*gdHy;

fluxHx =  gny.*gaEz +               alpha.*(gndotdH.*gnx-gdHx);
fluxHy = -gnx.*gaEz +               alpha.*(gndotdH.*gny-gdHy);
fluxEz = -gnx.*gdHy + gny.*gdHx   - alpha.*gdEz;
fluxHx = fluxHx/2;
fluxHy = fluxHy/2;
fluxEz = fluxEz/2;

rhsHx = rhsHx + (V*V'*gaussI')*((gauss.W.*fluxHx./gsqJ));
rhsHy = rhsHy + (V*V'*gaussI')*((gauss.W.*fluxHy./gsqJ));
rhsEz = rhsEz + (V*V'*gaussI')*((gauss.W.*fluxEz./gsqJ));

return;
