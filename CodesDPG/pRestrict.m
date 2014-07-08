% get interpolation matrix from N to Nr, acting on all elements 
function Irp = pRestrict(Norig,Nr)

Globals2D

[xo,yo] = Nodes2D(Norig); [ro, so] = xytors(xo,yo);

[xr,yr] = Nodes2D(Nr); [rr, sr] = xytors(xr,yr);
Vr = Vandermonde2D(Nr,rr,sr); 
Ir = Vandermonde2D(Nr,ro,so)/Vr;  % interp from Nr to Np points

Irp = kron(speye(K),Ir);
