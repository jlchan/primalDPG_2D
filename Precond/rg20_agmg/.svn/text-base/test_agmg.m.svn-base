%function x=test_agmg(A,Rhs,tol)
%matrixname = '';

%[A Rhs] = hdfread(matrixname);


% setup
tic
levels = agmg_setup(A);
toc
tol = 1e-4;
maxiter = 20;
tic
% solve
if (nargin<3)
    [x flag relres iter resvec] = agmg_solve(levels, Rhs);
else
    [x flag relres iter resvec] = agmg_solve(levels, Rhs, maxiter, tol);
end
toc


conv_rate = resvec(2:end)./resvec(1:end-1);
conv_rate = conv_rate(:);
conv_rate = [resvec(1)/norm(Rhs); conv_rate];

resvec = resvec/norm(Rhs);

    display(sprintf('iter  | conv.rate | rel res '));
    display(sprintf('-----------------------------'));
for i=1:length(resvec)
    display(sprintf('%4d  | %.4f | %.4e', i, conv_rate(i), resvec(i)));
end


display(sprintf('L2 residue = %e, relative residue = %e', norm(Rhs- ...
                                                  A*x), norm(Rhs-A* ...
                                                  x)/norm(Rhs)));
