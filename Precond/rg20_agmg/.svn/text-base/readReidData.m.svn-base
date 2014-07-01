function [A, f] = readReidData(filename)

A1_im = h5read(filename, '/A1_im');
A1_re = h5read(filename, '/A1_re');
I1 = h5read(filename, '/I1');
J1 = h5read(filename, '/J1');
f1 = h5read(filename, '/f1');


nnz = length(A1_im);
n = length(f1);

A_im = zeros(nnz,1);
A_re = zeros(nnz,1);
I = zeros(nnz,1);
J = zeros(nnz,1);
f = zeros(n, 1);


for i=1:nnz
    A_im(i) = A1_im(i);
    A_re(i) = A1_re(i);
    I(i) = I1(i);
    J(i) = J1(i);
end

for i=1:n
    f(i) = f1(i);
end

val = A_re + j*A_im;

A = spconvert([I J val]);