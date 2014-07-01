function [coarseA nullcoarseA P R] = coarsen(Amat, nullA);
% function [coarseA nullcoarseA P R] = coarsen(A, nullA)
% purpose : obtain coarse grid, interpolation and restriction operators
%           Assumes that A is given in matlab format


[m n] = size(Amat);

D = diag(Amat);
sign = (D>=0)-(D<0);

S = spdiags(-sign, 0, m, m);
SA = S*Amat;

% compute maximum off-diagonal
[dummy1 dummy2 maxOD] = find(max(SA(1:m,:))); maxOD = maxOD(:);

[I J V] = find(SA);

if(max(J) > n )
    display('error in cols of matrix');
end

if(max(I) > m )
    display('error in rows of matrix');
end



threshold = 0.5;

V = V - threshold*maxOD(I);

% keep the diagonals and strong neighbours
ids = find(V >= 0 | I == J);

C = sparse(I(ids), J(ids), V(ids));

C = mat2ell(C);

nnz_per_row = C.nnz_per_row;
nz = C.nz;
n = C.nrows;
ja = C.ja;
aa = C.aa;


% mis(2)
states = zeros(n,1);
r  = rand(n,1);

Ti = 1:1:n; Ti = Ti(:);
Tr = r;
Ts = states;

for i=1:n
    % dont count diagonal
    ids = find((ja(i,2:end) > 0));
    cols = ja(i,ids);
    Tr(cols) = Tr(cols)+1;
end

r = Tr;

Ti_hat = Ti; Tr_hat = Tr; Ts_hat = Ts;

done = false;

while(~done)

    for dist=1:2
        rmax = Tr;
        smax = Ts;
        imax = Ti;
        for jj = 1:nnz_per_row
            ids = find((ja(:,jj)>0));
            ids = ids(:);
            cols = ja(ids, jj);

            is_max = [];
            is_max = (smax(ids) > Ts(cols));

            is_max = is_max | ( (smax(ids) == Ts(cols)) & (rmax(ids) ...
                                                            > Tr(cols)) ...
                                 );
            is_max = is_max | ( (smax(ids) == Ts(cols)) & (rmax(ids) ...
                                                            == ...
                                                            Tr(cols)) ...
                                 & (imax(ids) > Ti(cols)));

            is_max = ~is_max;
            maxids = ids(find(is_max == 1));

            smax(maxids) = Ts(ja(maxids, jj));
            rmax(maxids) = Tr(ja(maxids, jj));
            imax(maxids) = Ti(ja(maxids, jj));
        end

        Ts = smax;
        Tr = rmax;
        Ti = imax;
    end

    % if states = 0 and imax = i
    ids = find((states == 0) & (Ti == (1:n)'));
    states(ids) = 1;

    % if states = 0 and smax = 1
    ids = find((states == 0) & (Ts == 1));
    states(ids) = -1;

    % reset tuples
    Ts = states;
    Tr = r;
    Ti = (1:n)';

    display(sprintf('number of undecided nodes = %d \n', ...
                    sum(states==0)));
    display(sprintf('number of removed nodes = %d \n', sum(states== ...
                                                      -1)));
    display(sprintf('number of MIS nodes = %d \n', sum(states==1)));
    % if at least one node is not decided
    done = (sum((states == 0)) == 0);
end

num_aggregates = sum(states == 1);

display(sprintf('number of aggregates = %d \n', num_aggregates));

% enumerate
agg = zeros(n,1);
agg = (states == 1);
agg = cumsum(agg);

if(agg(end) ~= num_aggregates)
    display('error in enumeration');
end

agg(states ~= 1) = 0;

Tr = r;
Ts = states;
Ti = (1:m)';
for dist=1:2
    rmax = Tr;
    smax = Ts;
    imax = Ti;
    for jj = 1:nnz_per_row
        ids = find((ja(:,jj)>0));
        ids = ids(:);
        cols = ja(ids, jj);
        is_max = [];
        is_max = (smax(ids) > Ts(cols));
        is_max = is_max | ( (smax(ids) == Ts(cols)) & (rmax(ids) ...
                                                        > Tr(cols)) ...
                             );
        is_max = is_max | ( (smax(ids) == Ts(cols)) & (rmax(ids) ...
                                                        == ...
                                                        Tr(cols)) ...
                             & (imax(ids) > Ti(cols)));

        is_max = ~is_max;
        maxids = ids(is_max);
        smax(maxids) = Ts(ja(maxids, jj));
        rmax(maxids) = Tr(ja(maxids, jj));
        imax(maxids) = Ti(ja(maxids, jj));
    end

    Ts = smax;
    Tr = rmax;
    Ti = imax;

    ids = find(agg == 0 & smax == 1 & agg(imax) ~= 0);
    agg(ids) = agg(imax(ids));
end

if(unique(sort(agg)) ~= (1:num_aggregates)')
    display('error in enumerating..');
end

P_ia = (1:n)';
P_ja = agg;

% P.aa = 0*agg + 1;
P_aa = nullA;

P = sparse(P_ia, P_ja, P_aa);

DP = sqrt(diag(P'*P));

nullcoarseA = DP;

P_aa = (P_aa)./nullcoarseA(P_ja);

P = sparse(P_ia, P_ja, P_aa);

R = P';

% Galerkin product
coarseA = R*Amat*P;

% condest(coarseA)
return;


function is_max = compare_tuple(s0, r0, i0, s1, r1, i1)
% returns 1 if (s0, r0, i0) >= (s1, r1, i1)
% returns 0 if (s0, r0, i0) < (s1, r1, i1)

if s0 > s1
    is_max = 1; return;
end

if s1 > s0
    is_max = 0; return;
end

if r0 > r1
    is_max = 1; return;
end

if r1 > r0
    is_max = 0; return;
end

if i0 > i1
    is_max = 1; return;
 end

if i1 > i0
    is_max = 0; return;
end

is_max = 1;

return;
