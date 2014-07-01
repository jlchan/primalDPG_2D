function rho = rhoDinvA(A)

% extract diagonal matrix
Diag = diag(A);
n = length(Diag);

Diag = spdiags(1./Diag, 0, n, n);

% form (D^-1) A
DinvA = Diag*A;

% estimate the largest eigen value from Arnoldi iteration
k = 5;
v0 = rand(n,1);
%v0 = ones(n,1);
[V H kact] = arnoldi(DinvA, v0, k, 0);
H = H(1:k,1:k);

rho = max(abs(eig(H)));

return