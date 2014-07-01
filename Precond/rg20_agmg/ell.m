function [s] = ell(m,n,nnz_per_row);

nz = nnz_per_row*m;

s = struct('nrows',m,'ncols',n,'nz', nz,'nnz_per_row', nnz_per_row, 'ja', zeros(m,nnz_per_row), ...
           'aa',zeros(m, nnz_per_row));

% nz = number of nonzero entries in the matrix
% m = number of rows
% n = number of columns
% ia = pointers to row starters
% ja = column indices
% aa = nonzero entries of matrix

return