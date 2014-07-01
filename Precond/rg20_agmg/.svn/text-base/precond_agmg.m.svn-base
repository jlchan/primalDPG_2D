function x = precond_agmg(levels, b);
% function x = precond_agmg(levels, b);
% Purpose : To obtain a solution of linear system via Kcycles
%
% Input  :
% levels : handle for AMG hierarchy
%     b  : rhs vector
% Output :
%     x  : solution vector



% rhs = b
levels{1}.rhs = b;

% start with zero initial guess
levels{1}.x = 0*b;

% call to recursive k-cycle at first level
%levels = vcycle(levels, 1);
levels = kcycle(levels, 1);
x = levels{1}.x;

%[x levels] = kcycle_notay(b, levels, 1);



% display('preconditioned...');
return;