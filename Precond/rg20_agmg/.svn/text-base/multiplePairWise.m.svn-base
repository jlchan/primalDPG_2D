function [coarseN coarseA AggIndex] = multiplePairWise(A, kTG, npass, tau)


[m n] = size(A);

[nc AggIndex] = initialPairWise(A, kTG);

P = sparse(m, nc);

cutoff = nnz(A)/tau;

for i=1:m
    if(AggIndex(i) > 0)
        P(i, AggIndex(i)) = 1;
    end
end
savedA = A;
A = P'*A*P;

for s=2:npass


    [nc AggIndex] = furtherPairWise(A, savedA, P, kTG);
    [m n] = size(A);

    P1 = sparse(m, nc);
    for i=1:m
        if(AggIndex(i) > 0)
            P1(i, AggIndex(i)) = 1;
        end
    end
    A = P1'*A*P1;

    P = P*P1;

    if(nnz(A) <= cutoff)
        break;
    end

end

coarseN = nc;
coarseA = A;

return;




function [nc AggIndex] = initialPairWise(A, kTG)

% initialize
[m n] = size(A);
AggIndex = -ones(m,1);

symA = 0.5*(A + A');

diagSymA = diag(symA);

S = sum(abs(symA), 2) - abs(diagSymA);

ids = find(diagSymA >= (kTG/(kTG-2))*S);

AggIndex(ids) = 0;

nc = 0;
%[J I symV] = find(symA');

while (sum(AggIndex == -1) > 0)

    for i=1:m
        if(AggIndex(i) == -1)

            % find nonzero column indices corresponding to row i
            [dummy J symV] = find(symA(i,:)); J = J(:); symV = symV(:);
            I = 0*J + i;

            numer = 2./(1./diagSymA(I) + 1./diagSymA(J) );
            denom = - symV + 1./(1./(diagSymA(I) - S(I)) + 1./(diagSymA(J) ...
                                                              - S(J)) );


            mu = numer./denom;
            mu = mu.*(mu > 0).*(diagSymA(I) - S(I) + diagSymA(J) - ...
                                S(J) >= 0);


            % remove already aggregated nodes
            ids = find(AggIndex(J) > -1); mu(ids) = 0;
            maxMu = 2*max(mu(mu >= 0));

            mu(mu == 0) = maxMu;

            [val id] = min(mu); j = J(id);

            nc = nc+1;

            if(val > 0 & val < maxMu)
                if(val <= kTG)
                    AggIndex(i) = nc; AggIndex(j) = nc;
                else
                    AggIndex(i) = nc;
                end
            else
                AggIndex(i) = nc;
            end


        end
    end

end
return



function [nc AggIndex] = furtherPairWise(A, savedA, P, kTG)

D = diag(savedA); m = length(savedA);
symA = 0.5*(savedA + savedA') - spdiags(D, 0, m, m);

test = P'*abs(symA)*P;

s = sum(test, 2);

symA = 0.5*(A+A');
D = diag(symA);
[n n] = size(A);
AggIndex = -ones(n, 1);

nc = 0;

while (sum(AggIndex == -1) > 0)

    for i=1:n

        if(AggIndex(i) == -1)

            % find nonzero column indices corresponding to row i
            [dummy J symV] = find(symA(i,:)); J = J(:); symV = symV(:);
            I = 0*J + i;


            numer = 2./(1./D(I) + 1./D(J) );
            denom = - symV + 1./(1./(D(I) - s(I)) + 1./(D(J) - s(J)) );


            mu = numer./denom;
            mu = mu.*(mu > 0).*(D(I) - s(I) + D(J) - s(J) >= 0);


            % remove already aggregated nodes
            ids = find(AggIndex(J) > -1); mu(ids) = 0;
            maxMu = 2*max(mu(mu >= 0));

            mu(mu == 0) = maxMu;

            [val id] = min(mu); j = J(id);

            nc = nc+1;
            if(val > 0 & val < maxMu)
                if(val <= kTG)
                    AggIndex(i) = nc; AggIndex(j) = nc;
                else
                    AggIndex(i) = nc;
                end
            else
                AggIndex(i) = nc;
            end


        end
    end

end