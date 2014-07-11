function [x flag relres iter resvec] = fpcg(A, b,tol, maxIter, Pre);

% grab size of matrix
[m n] = size(A);

% initial residue
nb = norm(b);

% initial guess
x0 = zeros(m, 1);

% initial residue
r0 = b - A*x0;

r = r0;
x = x0;

D = [];
AD = [];
mi = 1;

flag = 1;
resvec = zeros(maxIter, 1);

for i=1:maxIter

    % apply preconditioner (make sure that precond is symmetric)
    di = Pre(r);

    % without precondition
    % wi = r;

    rho = di'*r;

    if(i>1)
        % flexible conjugate gradient
        beta = (rho - di'*r_previous)/rho_previous;
        di = di + beta*d_previous;
    end

    Adi = A*di;

    alpha = (rho/(di'*Adi));

    % update solution
    x = x + alpha * di;

    r_previous = r;

    % update residue
    r = r - alpha * Adi;

    %    display(sprintf('residue = %e, actual residue = %e',
    %    norm(r), norm(b - A*x)));
    % display(sprintf('residue = %e', norm(r)));
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
relres = norm(b - A*x)/nb;
iter = maxIter;