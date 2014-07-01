function [zk levels] = kcycle_notay(rk, levels, k);

agmg_globals;

nLevels = length(levels);

if(k == nLevels) return; end

% zero out at the beginning
if(k > 1) levels{k}.x = 0*levels{k}.x; end

%display(sprintf('level = %d, res before pre-smooth = %e\n', k, norm(levels{k}.rhs)));

% pre smooth
zk1 = smooth(levels{k}.A, rk, 0*levels{k}.x, ...
             levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
             levels{k}.U);


% residual
rktilde = rk - levels{k}.A * zk1;

%display(sprintf('level = %d, res after pre-smooth = %e\n', k, norm(rktilde)));


% restrict the residual
rkp1 = levels{k}.R * rktilde;

if(k+1 < nLevels)

    % first krylov inner iteration
    [ckp1 levels] = kcycle_notay(rkp1, levels, k+1);

    vkp1 = levels{k+1}.A * ckp1;

    if strcmp(ktype, 'CG')
        rho1   = ckp1'*vkp1;
        alpha1 = ckp1'*rkp1;
    end

    if strcmp(ktype, 'GMRES')
        rho1   = vkp1'*vkp1;
        alpha1 = vkp1'*rkp1;
    end

    norm_rkp1 = norm(rkp1);
    rkp1_tilde = rkp1 - (alpha1/rho1)*vkp1;

    norm_rktilde_p = norm(rkp1_tilde);

    % verify this ..
    %    levels{k+1}.rhs = rkp1; % ??

    t = 0.2;
    if(norm_rktilde_p < t*norm_rkp1)
        % if first iteration is satisfactory update the solution
        % and get out of the loop
        xkp1_tilde = (alpha1/rho1)*ckp1;
    else
        % if first iteration is not satisfactory
        % do one more inner krylov iteration
        [dkp1 levels] = kcycle_notay(rkp1_tilde, levels, k+1);

        wkp1 = levels{k+1}.A * dkp1;

        if strcmp(ktype, 'CG')
            gam = dkp1'*vkp1;
            beta = dkp1'*wkp1;
            alpha2 = dkp1'*rkp1_tilde;
        end

        if strcmp(ktype, 'GMRES')
            gam = wkp1'*vkp1;
            beta = wkp1'*wkp1;
            alpha2 = wkp1'*rkp1_tilde;
        end

        % its tricky.. if rho1 is zero dont do anything its dangerous
        if(rho1 ~= 0)
            rho2 = beta - gam*gam/rho1;
            if(rho2 ~= 0)
                xkp1_tilde = (alpha1/rho1 - ((gam*alpha2)/(rho1*rho2)))*ckp1 ...
                    + (alpha2/rho2)*dkp1;
            end
        end
    end
else

    % for the coarsest grid use a direct solver
    xkp1_tilde = levels{k+1}.A\rkp1;
    %levels{k+1}.x = 0 * levels{k+1}.x;
end

zk2 = levels{k}.P * xkp1_tilde;

rkbar = rktilde - levels{k}.A * zk2;

% post smooth
zk3 = smooth(levels{k}.A, rkbar, 0*levels{k}.x, ...
             levels{k}.smoother, levels{k}.damping, levels{k}.L, ...
             levels{k}.U);

zk = zk1 + zk2 + zk3;

%norm(zk)

%display(sprintf('level = %d, res after post-smooth = %e\n', k, norm(rk - levels{k}.A*zk) ));
return;