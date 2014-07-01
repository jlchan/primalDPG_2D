function [x flag relres iter resvec] = agmg_solve(levels, b, maxIter, tol, x0);
% function [x flag relres iter resvec] = agmg_solve(levels, b);
% Purpose : Use Kcycles with GMRES/CG to obtain solution

agmg_globals;

display(sprintf('using %s iterative solver..', ktype));
display(sprintf('using %s smoother..', smoother_type));

% use zero as initial guess
if nargin < 5, x0 = 0*b;     end;
if nargin < 4, tol = 1e-6;   end;
if nargin < 3, maxIter = 40; end;


A = levels{1}.A;

% norm(b - A * x0)

% dont restart
% call for matlab gmres
% [x flag relres iter resvec] = gmres(A, b, maxIter, ...
%                                    tol, 1,
%                                    @(x)precond_agmg(levels,x));

% call for my gmres
if(strcmp(ktype, 'GMRES'))
    [x flag relres iter resvec] = my_gmres(A, b, maxIter, ...
                                           tol, 1, levels);
    return;
end
if(strcmp(ktype, 'CG'))
    [x flag relres iter resvec] = my_pcg(A, b, maxIter, ...
                                         tol, levels);
    return;
end