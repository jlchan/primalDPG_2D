function [x flag relres iter resvec] = my_pcg_nullA(A, b, maxIter, tol, levels);
% function [x flag relres iter resvec] = my_pcg(A, b, maxIter, tol, levels)
% Purpose : To solve a linear system iteratively with AMG
%           preconditioned generalized  minimal residual method
% Input   :
%     A   : matrix in matlab format
%     b   : right hand side vector
%  maxIter: maximum no of iteration
%    tol  : residual tolerance
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

nullA = levels{1}.nullA;


% initial guess
x0 = zeros(m, 1);

% initial residue
r0 = b - A*x0 - nullA*(nullA'*x0);

r = r0;
x = x0;

D = [];
AD = [];
mi = 1;

flag = 1;
resvec = zeros(maxIter, 1);

for i=1:maxIter

    % apply preconditioner (make sure that precond is symmetric)
    di = precond_agmg_nullA(levels, r);

    % without precondition
    % wi = r;

    rho = di'*r;

    if(i>1)
        % flexible conjugate gradient
        beta = (rho - di'*r_previous)/rho_previous;
        di = di + beta*d_previous;
    end

    Adi = A*di + nullA * (nullA'*di);


    alpha = (rho/(di'*Adi));

    % update solution
    x = x + alpha * di;

    r_previous = r;

    % update residue
    r = r - alpha * Adi;

    %    display(sprintf('residue = %e, actual residue = %e', norm(r), ...
    %                norm(b - A*x - nullA*(nullA'*x))));
    display(sprintf('residue = %e', norm(r)));
    resvec(i) = norm(r);

    d_previous = di;
    rho_previous = rho;

    % yikes.. dont do this
    % D = [D di];
    % AD = [AD Adi];

    if(resvec(i) < tol*nb)
        maxIter = i;
        flag = 0;
        break;
    end
end

%(D'*AD)./(D'*D)

resvec = resvec(1:maxIter);
relres = norm(b - A*x - nullA*(nullA'*x))/nb;
iter = maxIter;