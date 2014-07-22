function [IM] = InterpMatrix2D(rout, sout)

% function [IM] = InterpMatrix2D(rout, sout)
% purpose: compute local elemental interpolation matrix

Globals2D;
 
% compute Vandermonde at (rout,sout)
Vout = Vandermonde2D(N, rout, sout);

% build interpolation matrix (Interp = Vout*inv(V) => Interp*V*V' = Vout*V' => Interp = Vout*V'/(V*V')
IM = Vout/V;

return

  
  
