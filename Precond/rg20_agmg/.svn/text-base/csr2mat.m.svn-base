function Amat = csr2mat(A_csr)
% function Amat = csr2mat(A_csr)
% purpose : converts a matrix stored in csr format to matlab format

n = A_csr.nrows;
m = A_csr.ncols;
ia = A_csr.ia;
js = A_csr.ja;
aa = A_csr.aa;

is = 0*js;

for i=1:n
    Jstart = ia(i);
    Jend = ia(i+1)-1;
    for jj=Jstart:Jend
        is(jj) = i;
    end
end


Amat = sparse(is, js, aa, n, m);