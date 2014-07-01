function [coarseA nullcoarseA P R] = coarsen(Amat, nullA);
% function [coarseA nullcoarseA P R] = coarsen(A, nullA)
% purpose : obtain coarse grid, interpolation and restriction operators
%           Assumes that A is given in matlab format


A = mat2csr_opt(Amat);

nz = A.nz;
n = A.nrows;
ia = A.ia;
ja = A.ja;
aa = A.aa;

threshold = 0.5;

is_strong = zeros(nz,1);

%sign = (aa(ia(1:n)) > 0) - aa(ia(1:n)<0);
%sign = sign(:);

% TODO: vectorize
for i=1:n
    Jstart = ia(i);
    Jend = ia(i+1)-1;

    if(ja(Jstart) ~= i)
        display('error in csr format');
        return;
    end

    max_OD = 0;
    sign = (aa(Jstart) > 0) - (aa(Jstart) < 0);

    for jj=Jstart+1:Jend
        if(-sign*aa(jj) > max_OD)
            max_OD = -sign*aa(jj);
        end
    end

    is_strong(Jstart) = 1;
    for jj=Jstart+1:Jend
        is_strong(jj) = (-sign*aa(jj) >= threshold*max_OD);
    end
end

% mis(2)
states = zeros(n,1);
r  = rand(n,1);

Ti = 1:1:n; Ti = Ti(:);
Tr = r;
Ts = states;


for i=1:nz
    Tr(ja(i)) = Tr(ja(i)) + 1;
end

Ti_hat = Ti; Tr_hat = Tr; Ts_hat = Ts;

done = false;

while(~done)

    for dist=1:2
        for i=1:n
            Jstart = ia(i);
            Jend = ia(i+1)-1;
            smax = Ts(i); rmax = Tr(i); imax = Ti(i);

            for jj=Jstart+1:Jend
                if is_strong(jj)
                    c = ja(jj);
                    sc = Ts(c); rc = Tr(c); ic = Ti(c);
                    if(compare_tuple(smax, rmax, imax, sc, rc, ic) == 0)
                        smax = sc; rmax = rc; imax = ic;
                    end
                end
            end

            Ts_hat(i) = smax;
            Tr_hat(i) = rmax;
            Ti_hat(i) = imax;
        end

        Ts = Ts_hat;
        Tr = Tr_hat;
        Ti = Ti_hat;
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

% TODO : vectorize using cumsum
% enumerate
agg = zeros(n,1);
aggIndex = 1;
for i=1:n
    if states(i) == 1
        agg(i) = aggIndex;
        aggIndex = aggIndex + 1;
    end
end


% mis to aggregates
for dist=1:2
    for i=1:n
        Jstart = ia(i);
        Jend = ia(i+1)-1;
        smax = Ts(i); rmax = Tr(i); imax = Ti(i);

        for jj=Jstart+1:Jend
            if is_strong(jj)
                c = ja(jj);
                sc = Ts(c); rc = Tr(c); ic = Ti(c);
                if(compare_tuple(smax, rmax, imax, sc, rc, ic) == 0)
                    smax = sc; rmax = rc; imax = ic;
                end
            end
        end

        Ts_hat(i) = smax;
        Tr_hat(i) = rmax;
        Ti_hat(i) = imax;

        if(agg(i) == 0 && smax == 1 && agg(imax) ~= 0)
            agg(i) = agg(imax);
        end
    end

    Ts = Ts_hat;
    Tr = Tr_hat;
    Ti = Ti_hat;

    %    agg = agg(Ti);
end

if(unique(sort(agg)) ~= (1:num_aggregates)')
    display('error in enumerating..');
end

P = csr(n, num_aggregates, n);
P.nrows = n;
P.ncols = num_aggregates;
P.nz = n;
P.ia = (1:n+1)';
P.ja = agg;

% P.aa = 0*agg + 1;
P.aa = nullA;

Pmat = csr2mat(P);

DP = sqrt(diag(Pmat'*Pmat));

nullcoarseA = DP;



P.aa = (P.aa)./nullcoarseA(P.ja);

P = csr2mat(P);

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
