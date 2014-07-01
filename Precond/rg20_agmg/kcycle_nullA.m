function [levels] = kcycle_nullA(levels, k);

agmg_globals;

nLevels = length(levels);

if(k == nLevels) return; end

nullAk = levels{k}.nullA;

nullAkp1 = levels{k+1}.nullA;

% zero out at the beginning
if(k > 1) levels{k}.x = 0*levels{k}.x; end

%display(sprintf('level = %d, res before pre-smooth = %e\n', k, norm(levels{k}.rhs)));

% pre smooth
levels{k}.x = smooth_nullA(levels{k}.A, nullAk, levels{k}.rhs, levels{k}.x, ...
                     levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
                     levels{k}.U);


% residual
res = levels{k}.rhs - levels{k}.A * levels{k}.x - nullAk*(nullAk'*levels{k}.x);


%display(sprintf('level = %d, res after pre-smooth = %e\n', k, norm(res)));


% restrict the residual
levels{k+1}.rhs = levels{k}.R * res;

if(k+1 < nLevels)

    % first krylov inner iteration
    levels = kcycle_nullA(levels, k+1);

    ckp1 = levels{k+1}.x;
    vkp1 = levels{k+1}.A * ckp1 +  nullAkp1 * (nullAkp1'*ckp1);

    rkp1 = levels{k+1}.rhs;

    if strcmp(ktype, 'CG')
        rho1   = ckp1'*vkp1;
        alpha1 = ckp1'*rkp1;
    end

    if strcmp(ktype, 'GMRES')
        rho1   = vkp1'*vkp1;
        alpha1 = vkp1'*rkp1;
    end

    norm_rkp1 = norm(rkp1);
    rkp1 = rkp1 - (alpha1/rho1)*vkp1;

    norm_rktilde_p = norm(rkp1);

    % verify this ..
    levels{k+1}.rhs = rkp1; % ??

    t = 0.2;
    if(norm_rktilde_p < t*norm_rkp1)
        % if first iteration is satisfactory update the solution
        % and get out of the loop
        levels{k+1}.x = (alpha1/rho1)*ckp1;
    else
        % if first iteration is not satisfactory
        % do one more inner krylov iteration
        levels = kcycle_nullA(levels, k+1);
        dkp1 = levels{k+1}.x;

        wkp1 = levels{k+1}.A * dkp1 + nullAkp1*(nullAkp1'*dkp1);

        if strcmp(ktype, 'CG')
            gam = dkp1'*vkp1;
            beta = dkp1'*wkp1;
            alpha2 = dkp1'*rkp1;
        end

        if strcmp(ktype, 'GMRES')
            gam = wkp1'*vkp1;
            beta = wkp1'*wkp1;
            alpha2 = wkp1'*rkp1;
        end

        % its tricky.. if rho1 is zero dont do anything its dangerous
        if(rho1 ~= 0)
            rho2 = beta - gam*gam/rho1;
            if(rho2 ~= 0)
                levels{k+1}.x = (alpha1/rho1 - ((gam*alpha2)/(rho1*rho2)))*ckp1 ...
                    + (alpha2/rho2)*dkp1;
            end
        end
    end
else

    M = levels{k+1}.A + nullAkp1*nullAkp1';
    % for the coarsest grid use a direct solver
    levels{k+1}.x = M\levels{k+1}.rhs;
    %levels{k+1}.x = 0 * levels{k+1}.x;
end

% prolongate
levels{k}.x = levels{k}.x + levels{k}.P * levels{k+1}.x;

%display(sprintf('level = %d, res before post-smooth= %e\n', k, norm(levels{k}.rhs - levels{k}.A*levels{k}.x) ));


% post smooth
levels{k}.x = smooth_nullA(levels{k}.A, nullAk, levels{k}.rhs, levels{k}.x, ...
                           levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
                           levels{k}.U);

%display(sprintf('level = %d, res after post-smooth = %e\n', k, norm(levels{k}.rhs - levels{k}.A*levels{k}.x) ));
return;