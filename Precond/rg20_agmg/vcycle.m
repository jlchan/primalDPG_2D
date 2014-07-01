function [levels] = vcycle(levels, k);

nLevels = length(levels);

if(k == nLevels) return; end

% zero out at the beginning
if(k > 1) levels{k}.x = 0*levels{k}.x; end

% display(sprintf('level = %d, res before pre-smooth = %e\n', k,
% norm(levels{k}.rhs)));

% pre smooth
levels{k}.x = smooth(levels{k}.A, levels{k}.rhs, levels{k}.x, ...
                     levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
                     levels{k}.U);


% residual
res = levels{k}.rhs - levels{k}.A * levels{k}.x;

% display(sprintf('level = %d, res after pre-smooth = %e\n', k, norm(res)));

% restrict the residual
levels{k+1}.rhs = levels{k}.R * res;

if(k+1 < nLevels)
    levels = vcycle(levels, k+1);
else

    % for the coarsest grid use a direct solver
    levels{k+1}.x = levels{k+1}.A\levels{k+1}.rhs;
    %    levels{k+1}.x = 0 * levels{k+1}.x;
end

% prolongate
levels{k}.x = levels{k}.x + levels{k}.P * levels{k+1}.x;

% display(sprintf('level = %d, res before post-smooth = %e\n', k, norm(levels{k}.rhs - levels{k}.A*levels{k}.x) ...
%              ));


% post smooth
levels{k}.x = smooth(levels{k}.A, levels{k}.rhs, levels{k}.x, ...
                     levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
                     levels{k}.U);
return;