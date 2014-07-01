function [x flag relres iter resvec] = my_gmres(A, b, maxIter, tol, outIter, levels);
% function [x flag relres iter resvec] = my_gmres(A, b, maxIter, tol, ...
%                                           outerIter, levels)
% Purpose : To solve a linear system iteratively with AMG
%           preconditioned generalized  minimal residual method
% Input   :
%     A   : matrix in matlab format
%     b   : right hand side vector
%  maxIter: maximum no of iteration
%    tol  : residual tolerance
% outIter : outer iterations for restarted gmres (ignored for this code)
% levels  : AMG hierarchy
%
% Output  :
%      x  : solution vector
%   flag  : is zero if the solution is converged to given tol
%  relres : final relative residual L2 norm( b - A*x)/norm(b)
%    iter : number of iterations used
%  resvec : history of the residual for each iteration



% grab size of matrix
[m n] = size(A);

% initial residue
nb = norm(b);

% initial guess
u0 = zeros(m, 1);

% initial residue
r0 = b - A*u0;

r = r0;

H = zeros(maxIter, maxIter);
alphas = zeros(maxIter, 1);
resvec = zeros(maxIter, 1);
flag = 1;

%
C = [];
Z = [];

for j=1:maxIter
    % apply preconditioner
    z = precond_agmg(levels, r);

    % matvec
    cj = A*z;

    c = cj;

    % orthogonalize c w.r.t to previous vectors
    for i=1:j-1

        ci = C(:,i);

        gij = ci'*cj; H(i,j) = gij;

        c = c - gij*ci;
    end

    % normalize 'c'
    gjj = norm(c); H(j,j) = gjj;

    if abs(gjj) < tol*nb
        flag = 3;
        display('gjj is too small.. stopping');
        maxIter = j-1;
        break;
    end

    c = c/gjj;

    % update residual
    alpha_j  = r'*c; r = r - alpha_j*c;

    alphas(j) = alpha_j;

    resvec(j) = norm(r);

    % yikes.. dont do this
    C = [C c];
    Z = [Z z];

    if norm(r) < tol*nb
        flag = 0;
        maxIter = j;
        break;
    end

end


H = H(1:maxIter, 1:maxIter);
alphas = alphas(1:maxIter);
resvec = resvec(1:maxIter);

% get the solution
x = Z*(H\alphas);

% compute rel res
relres = norm(b - A*x)/nb;

iter = maxIter;