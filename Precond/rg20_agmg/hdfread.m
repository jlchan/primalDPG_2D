function [A Rhs] = hdfread(matrixname)

h5name = matrixname;
%h5name = sprintf('%s.h5', matrixname);

num_rows = h5read(h5name, '/num_rows');
num_cols = h5read(h5name, '/num_cols');

rowptrs = h5read(h5name, '/rowptrs');
cols = h5read(h5name, '/cols');
coefs = h5read(h5name, '/coefs');

rhs = h5read(h5name, '/rhs');


nz = length(coefs);

% since they are in c-format
if((length(rowptrs) ~= num_rows+1) || (rowptrs(end) ~= nz))
    display('error in rowptrs');
end


rowptrs = rowptrs + 1;

Rows = zeros(nz,1);
Cols = zeros(nz,1);
Coefs = zeros(nz,1);
Rhs = zeros(num_cols,1);

for i=1:num_rows
    Rows(rowptrs(i):rowptrs(i+1)-1) = i;
    Rhs(i) = rhs(i);
end

for i=1:nz
    Cols(i) = cols(i)+1;
    Coefs(i) = coefs(i);
end

A = sparse(Rows, Cols, Coefs);