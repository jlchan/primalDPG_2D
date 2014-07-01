function [coarseA nullcoarseA P R] = sag_coarsen(Amat, nullA);
% function [coarseA nullcoarseA P R] = coarsen(A, nullA)
% purpose : obtain coarse grid, interpolation and restriction operators
%           Assumes that A is given in matlab format


agmg_globals;

[m n] = size(Amat);

D = diag(Amat);
D = abs(D).^(-0.5);
D = spdiags(D, 0, m, m);

S = D*abs(Amat)*D;

C = (S > 0.08);

[CI CJ CV] = find(C);

% mis(2)
% state = 0 removed
%       = 1 undecided
%       = 2 MIS node
states = ones(m,1);
r  = rand(m,1);

rand_max = ceil(max(r))+1;


done = false;

while(~done)

    my_key = states*rand_max + r;
    max_key = my_key;
    for dist=1:2

        %        max_key = Ts*rand_max + r;
        key = sparse(CI, CJ, max_key(CJ));

        % find max key in each row
        [d1 d2 max_key(:)] = find(max(key(1:m,:), [], 2));

        max_key = max_key(:);

    end

    smax = floor(max_key/rand_max);

    % MIS nodes
    ids = find((states == 1) & (max_key == my_key));
    states(ids) = 2;

    % removed nodes
    ids = find(states == 1 & smax == 2);
    states(ids) = 0;

    % display(sprintf('number of undecided nodes = %d \n', ...
    %                 sum(states==1)));
    % display(sprintf('number of removed nodes = %d \n', sum(states== ...
    %                                                   0)));
    % display(sprintf('number of MIS nodes = %d \n', sum(states==2)));
    % if at least one node is not decided
    done = (sum((states == 1)) == 0);
end

states = states/2;

num_aggregates = sum(states == 1);

% display(sprintf('number of aggregates = %d \n', num_aggregates));

% enumerate
agg = zeros(n,1);
agg = (states == 1);
agg = cumsum(agg);

if(agg(end) ~= num_aggregates)
    display('error in enumeration');
end

agg(states ~= 1) = 0;

states = 2*states;
% aggregate non MIS nodes to mis nodes
my_key = states*rand_max + r;
max_key = my_key;
for dist = 1:2

    key = sparse(CI, CJ, max_key(CJ));

    % find max_key for each row
    [dummy1 dummy2 max_key(:)] = find(max(key(1:m,:), [], 2));
    [key_i key_j key_v] = find(key);

    max_key = max_key(:);

    ids = find(key_v == max_key(key_i));

    % extract column indices corresponding to max neigh for each
    % row
    max_cols = zeros(m,1);
    max_cols(key_i(ids)) = key_j(ids);

    if(sort(unique(key_i(ids))) ~= (1:m)')
        display('inappropriate');
    end

    smax = floor(max_key/rand_max);
    imax = max_cols;

    % size(find(agg==0))
    % size(find(smax == 2))
    % size(find(agg(imax)~=0))
    ids = find(agg == 0 & smax == 2 & agg(imax) ~= 0);

    agg(ids) = agg(imax(ids));
end

if(sum(agg == 0) > 0)
    display('danger..');
    display(sprintf('no of unaggregated nodes = %d\n', sum(agg==0)));
end

if(unique(sort(agg)) ~= (1:num_aggregates)')
    display('error in enumerating..');
end

P_ia = (1:m)';
P_ja = agg;

% P.aa = 0*agg + 1;
P_aa = nullA;

P = sparse(P_ia, P_ja, P_aa);

% compute column norms of P
DP = sqrt(diag(P'*P));

nullcoarseA = DP;

P_aa = (P_aa)./nullcoarseA(P_ja);

P = sparse(P_ia, P_ja, P_aa);


D = diag(Amat);
D = 1./D;
D = spdiags(D, 0, m , m);
DinvA = D*Amat;

rho_DinvA = rhoDinvA(Amat);

omega = (4.0/3.0)/rho_DinvA;

P = (spdiags(ones(m,1), 0, m, m) - omega*DinvA)*P;
R = P';


% norm(R*P - speye(num_aggregates), 1)

% Galerkin product
coarseA = R*Amat*P;

% condest(coarseA)
return;