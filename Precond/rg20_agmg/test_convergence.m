function rate = test_convergence(A, levels)
[m,n] = size(A);

rhs = rand(m,1);

rate = eigs(@(x)amg_error(x,levels), m, 1, 'LM');

display(sprintf('convergence rate of agmg = %e', rate));


% (I-M*A)*X
function Err = amg_error(x, levels)

A = levels{1}.A;

b = A*x;

y = precond_agmg(levels, b);

Err = x - y;
return;