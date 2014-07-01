function y = matvec_csr(A, x);
% function y = matvec_csr(A, x);
% Purpose : compute the product of matrix in CSR format with a
% vector

m = A.nrows;
n = A.ncols;

if(n ~= length(x))
    display('mismatch of dimensions..');
    return;
end

y  = zeros(m,1);
ia = A.ia;
ja = A.ja;
aa = A.aa;

for i = 1:n
    Jstart = ia(i);
    Jend   = ia(i+1)-1;
    y(i)   = (aa(Jstart:Jend))'*x(ja(Jstart:Jend));
end