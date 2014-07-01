function [s] = csr(m,n,nz);

s = struct('nrows',m,'ncols',n,'nz', nz,'ia',[1:m+1],'ja',[1:nz],'aa',rand(1,nz));
% nz = number of nonzero entries in the matrix
% m = number of rows
% n = number of columns
% ia = pointers to row starters
% ja = column indices
% aa = nonzero entries of matrix

return