function A_csr = mat2csr_opt(A)
% MAT2CSR   Converts Matlab-Matrix into compressed sparse rows
% format
%
%   [rst, rl, cn, val] = MAT2CSR(A) converts the Matlab-Matrix A
%   into the
%   compressed sparse rows format rst, rl, cn, val.
%
%   Note: This function makes use of the fact that Matlab stores
%   matrices
%   internally in a column oriented way.
%
%   See also find, sparse, accumarray, cumsum
%
%   David Fritzsche, 26 May 2008

% Note that the column oriented structure of A' is the row oriented
% structure of A.
% The column numbers (here stored in I) returned by find are
% growing
% monotone, i.e., in cn and val comes first the data for the first
% column
% of A' (first row of A), followed by the data for the second
% column of A'
% (second row of A), and so on.
[cn,I,val] = find(A');

% We determine the row lengths by using Matlabs powerful accumarray
% function.
rl = accumarray( I, 1, [size(A,1), 1] );

% The row-starts are just a cumulative sum of the row length. The
% vector of
% ones just fixes the indexing.
rst = cumsum([1; rl(1:end)]);

[m n] = size(A);

if(m == n)
    for i=1:m

        J = cn(rst(i):rst(i+1)-1);
        AA = val(rst(i):rst(i+1)-1);
        temp = abs(J(:)-i);
        [vals ids] = sort(temp);
        J = J(ids);
        AA = AA(ids);

        cn(rst(i):rst(i+1)-1) = J;
        val(rst(i):rst(i+1)-1) = AA;
    end
end

A_csr = csr(m,n,nnz(A));
A_csr.ia = rst;
A_csr.ja = cn;
A_csr.aa = val;
