% gets the boundary matrix Mb and volume to boundary matrix Eb - 
% the integral of u*v over the boundary is then v'*R*Eb'*Mb*Eb*R'*u.
%
% input: bmask = # face nodes, tells which face nodes to not put an
% integral on. if not specified, defaults to entire boundary

function [Mb Eb] = getBoundaryMatrix(bmask)

Globals2D

NfpB = length(vmapB);
nbfaces = NfpB/Nfp; 
Eb = sparse(1:NfpB,vmapB,1,NfpB,Np*K); % maps volume nodes to boundary nodes

% if not specified, assume we want integral over the entire boundary
if nargin<1
    bmask = ones(NfpB,1);
end

% make edge integration matrix on non-boundary edges
Nf = Nfp-1; r1D = JacobiGL(0,0,Nf); V1D = Vandermonde1D(Nf,r1D);
M1D = inv(V1D*V1D'); % 1D stiffness matrix from GLL nodes for faces

% note: may not need bfaces most of the time - setting BCs automatically
% gets rid of BC entries anyways, so we can possibly precompute Mb.
Mbface = spdiag(sJ(mapB))*kron(speye(nbfaces),M1D);
Mb = Mbface*spdiag(bmask(:)); 

return

