% returns Rr such that Rr*Ak*Rr' gives back a CG matrix of order Nr from the
% block diagonal matrix Ak (of order N >= Nr). 
% other return values: boundary indices (into reduced system) and nodal coordinates.
% - vmapBT, xrestrict, yrestrict 
function [Rp Irp vmapBT xrestrictB yrestrictB xrestrict yrestrict] = pRestrictCG(Norig,Nr)

Globals2D

saveFlag = 1;
a = backupGlobals(saveFlag);

% ugly way to do this part... fix - do a smarter restriction? or do no
% restriction at all, and just do elem computation locally. next time.
Nold = N; N = Nr; %save old N 
StartUp2D; % FaceStartUp2D; 

Rp = getCGRestriction();

[i j] = find(round(Rp));
% give back new boundary nodes
vmapBT = i(vmapB);
xrestrictB = x(vmapB);
yrestrictB = y(vmapB);

xrestrict = x;
yrestrict = y;

% reset after getting Rp
N = Nold;
% StartUp2D; 

saveFlag = 0;
backupGlobals(saveFlag,a);

Irp = pRestrict(Norig,Nr);
% Rr = Rp*Irp';
