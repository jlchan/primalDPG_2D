function [coarseA nullcoarseA P R] = coarsen(Amat, nullA);
% function [coarseA nullcoarseA P R] = coarsen(A, nullA)
% purpose : obtain coarse grid, interpolation and restriction operators
%           Assumes that A is given in matlab format


agmg_globals;

[m n] = size(Amat);

D = diag(Amat);
sign = (D>=0)-(D<0);

S = spdiags(-sign, 0, m, m);
SA = S*Amat;

% TODO : delete this work around .. and just use max

% a work around for matlab
% if a row has only one nonzero stick in one more zero
% so that max returns at least one value per row
[c I v] = find(SA');
%rl = accumarray(I, 1, [m 1]);
%ids1 = find(rl == 1);
%ids2 = find(rl > 1);
maxOD = zeros(m,1);


% compute maximum off-diagonal
%[dummy1 dummy2 maxOD(ids2)] = find(max(SA(1:m,:), [], 2 ) ); % mD
%= mD(:);
[dummy1 dummy2 mD] = find(max(SA(1:m,:), [], 2 ) ); % mD = mD(:);

maxOD(dummy1) = mD(:);

[I J V] = find(SA);

if(max(J) > n )
    display('error in cols of matrix');
end

if(max(I) > m )
    display('error in rows of matrix');
end


V = V - threshold*maxOD(I);

% keep the diagonals and strong neighbours
ids = find( (V > 0 & maxOD(I) > 0)| I == J );

% graph of strong connections
C = sparse(I(ids), J(ids), ones(length(ids), 1));

[CI CJ CV] = find(C);

% mis(2)
% state = 0 removed
%       = 1 undecided
%       = 2 MIS node
states = ones(m,1);
r  = rand(m,1);

%%%%%%r = (1:m)';

% r(j) += #{C_ij != 0}
% i.e. segmented reduction based on J
r = r + accumarray(CJ, 1, [m, 1]);

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

R = P';


% norm(R*P - speye(num_aggregates), 1)

% Galerkin product
coarseA = R*Amat*P;

% condest(coarseA)
return;