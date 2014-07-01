function levels = agmg_setup_nullA(A, nullA)
% function levels = agmg_setup(A)
% Purpose : construct agmg levels
%
% Input  :
%   A    : matrix in matlab format
%
% Output :
% levels : matlab object for stroring hierarchy of amg
%          levels{k} has info of level k
%          each level has A, R, P, smoother infos
%

agmg_globals;

% grab the size of matrix
[m n] = size(A);


if m ~= n
    display('matrix is not square');
    return;
end

if nargin < 2, nullA = ones(n,1); end

nullA = nullA/norm(nullA);

k = 1;

smoother = smoother_type;

% size of the coarse grid solver
nL = coarse_size;

while  (size(A, 1) >  nL)
    % construct coarse level matrix, interpolation operator,
    % and restriction operator


    % unsmooth aggregation
    [coarseA nullcoarseA P R] = coarsen(A, nullA);

    % smooth aggregation
    %[coarseA nullcoarseA P R] = sag_coarsen(A, nullA);

    [m n] = size(A);

    damping = 1;

    if strcmp (smoother, 'J')
        L = [];
        U = [];
    end
    if strcmp (smoother, 'symDJ')
        damping = symDJfactor_nullA(A, nullA);
        L = [];
        U = [];
    end
    if strcmp (smoother, 'DJ')
        damping = DJfactor_nullA(A, nullA);
        display(sprintf('Estimated spectral radius = %e', 4.0/(3.0*damping)));
        L = [];
        U = [];
    end

    if strcmp (smoother, 'GS')
        damping = 1;
        L = tril(A+nullA*nullA',0);
        U = triu(A+nullA*nullA',1);
    end

    if strcmp (smoother, 'symGS')
        damping = 1;
        L = tril(A+nullA*nullA',0);
        U = triu(A+nullA*nullA',0);
    end

    if strcmp (smoother, 'Cheb')
        damping = DJfactor_nullA(A, nullA);
        % use spectral radius instead
        damping = 4.0/(3.0*damping);
        L = [];
        U = [];
    end


    % copy the matrices to levels{k}
    levels{k} = struct('A', A, 'P', P, 'R', R, 'nullA', nullA, 'L', L, 'U',U, 'smoother', smoother, ...
                       'damping', damping, 'x', zeros(n,1), ...
                       'rhs', zeros(n,1), 'res', zeros(n,1));


    A = coarseA;
    nullA = nullcoarseA;

    nullA = nullA / norm(nullA);
    k = k+1;

    if(size(levels{k-1}.A, 1) < 2*size(A,1))
        break;
    end
end

[m n] = size(A);

L = [];
U = [];
% for final level
levels{k} = struct('A', A, 'P', zeros(1,1), 'R', zeros(1,1), 'nullA', ...
                   nullA, 'L', L, 'U', U, 'smoother', smoother, 'damping', 1, 'x', zeros(n,1), 'rhs', zeros(n,1), 'res', zeros(n,1));

numLevels = k;


% print the hierarchy
    display(sprintf('level |  rows   |  cols   |   nnz    |nnz/row| maxnnz |'));
    display(sprintf('------|---------|---------|----------|-------|--------|'));
for i=1:k

    [nrows ncols] = size(levels{i}.A);
    nz = nnz(levels{i}.A);

    [c I v] = find((levels{i}.A)');
    rl = accumarray(I, 1, [nrows 1]);

    display(sprintf('%4d  | %7d | %7d | %8d | %4.2f | %6d |', i, nrows, ncols, ...
                    nz, nz/nrows, max(rl)));
end


return



function damping = DJfactor(A)

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

% estimate the damping factor using
damping = (4.0/3.0)/rho;

return


function damping = DJfactor_nullA(A, nullA)

% extract diagonal matrix
n = length(nullA);

% estimate the largest eigen value from Arnoldi iteration
k = 5;
v0 = rand(n,1);
%v0 = ones(n,1);
[V H kact] = arnoldi_nullA(A, nullA, v0, k, 0);
H = H(1:k,1:k);

rho = max(abs(eig(H)));

% estimate the damping factor using
damping = (4.0/3.0)/rho;

return


function damping = symDJfactor(A)

% extract diagonal matrix
Diag = diag(A);
n = length(Diag);

D = spdiags(1./sqrt(abs(Diag)), 0, n, n);

% form (D^{-1/2}) A (D^{-1/2})
DAD = D*A*D;

k = 5;
[V H kact] = arnoldi(DAD, rand(n,1), k, 0);
H = H(1:k,1:k);


rho = max(abs(eig(H)));

% estimate the damping factor using
damping = (4.0/3.0)/rho;

return