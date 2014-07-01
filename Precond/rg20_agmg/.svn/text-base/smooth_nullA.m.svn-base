function x = smooth_nullA(A, nullA, b, x, type, dfactor, L, U);
% function x = smooth(A, b, x, dfactor);
% Purpose : Applies a smoother

if nargin < 5, dfactor = 1.0, end;
if nargin < 4, type = 'GS', end;

n = size(A, 1);

agmg_globals;

if nargin < 3, x = zeros(n, 1), end;

% if Jacobi smoother is being called
% x = x + D^{-1} (b - Ax)
if strcmp(type, 'J')
    D = diag(A);
    D = D + nullA.*nullA;
    r = b - A*x - nullA*(nullA'*x);
    x = x + (r./D);
    return;
end

% if damped Jacobi smoother is being called
% x = x + omega * D^{-1} (b - Ax)
if strcmp(type, 'DJ')
    D = diag(A);
    D = D + nullA.*nullA;
    r = b - A*x - nullA*(nullA'*x);
    x = x + dfactor * (r./D);
    return;
end


% if symmetric damped Jacobi smoother is being called
% x = x + omega * (D^{-1} b - D^{-1/2}AD^{-1/2} x)
if strcmp(type, 'symDJ')
    D = diag(A);
    D = D + nullA.*nullA;
    sqrtD = sqrt(abs(D));

    r = b - sqrtD.*(A*(x./sqrtD)) - sqrtD.*(nullA*(nullA'*(x./sqrtD)));
    x = x + dfactor * (r./D);
    return;
end

% if Gauss Seidel smoother is being called
% Lx = b - Ux
if strcmp(type, 'GS')
    % lower triangle
    %    L = tril(A, 0);
    % strict upper triangle
    %    U = triu(A, 1);

    r = b - U*x;

    x = L\r;

    return;
end


% if symmetric Gauss Seidel smoother is being called
if strcmp(type, 'symGS')

    % M = low(A)*invD*upp(A),
    %    L = tril(A, 0);
    %    U = triu(A, 0);
    D = diag(A) + nullA.*nullA;

    r = b - A*x - nullA*(nullA'*x);

    % x = inv(M) * r = inv(U)*D*inv(L)*r
    x = x + U\(D.*(L\r));
    return;
end


if strcmp(type, 'Cheb')

    D = diag(A) + nullA.*nullA;

    invD = 1./D;

    lambda_max = dfactor;

    lambda_min = lambda_max/10.;

    rho = (lambda_max - lambda_min)/(lambda_max + lambda_min);

    alpha = 0.25*rho*rho;

    gamma = 1;

    lambda_avg = 0.5*(lambda_min+lambda_max);

    invD = invD*(1/lambda_avg);

    xnm1 = 0*x;
    rnm1 = b - A*x - nullA*(nullA'*x);

    b = rnm1;

    xn = invD.*rnm1;

    gamma_n = 1;
    beta_n = 2;

    for i=1:chebyshev_order

        rn = b - A*xn - nullA*(nullA'*xn);

        gamma_np1 = alpha*beta_n/(1 - alpha*beta_n);

        beta_np1 = 1 + gamma_np1;

        xnp1 = beta_np1*(xn + invD.*rn) - gamma_np1*xnm1;

        xnm1 = xn;
        xn = xnp1;

        beta_n = beta_np1;
        gamma_n = gamma_np1;
    end

    x = x+xn;
    return;
end
return;