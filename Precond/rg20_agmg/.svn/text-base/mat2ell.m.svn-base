function A_ell = mat2ell(A)

[cn,I,val] = find(A');

% We determine the row lengths by using Matlabs powerful accumarray
% function.
rl = accumarray( I, 1, [size(A,1), 1] );

maxrl = max(rl(:));

% The row-starts are just a cumulative sum of the row length. The
% vector of
% ones just fixes the indexing.
rst = cumsum([1; rl(1:end)]);

[m n] = size(A);

A_ell = ell(m,n,maxrl);

display(sprintf('max row size = %d, avg row size = %d \n', maxrl, ...
                sum(rl)/m));

for i=1:m

    J = cn(rst(i):rst(i+1)-1);
    AA = val(rst(i):rst(i+1)-1);

    if(n==m)
        temp = -1e9*(J(:)==i) + J(:).*(J(:)~=i);
        [vals ids] = sort(temp);
        J = J(ids);
        AA = AA(ids);
    end

    %    display(sprintf('i = %d, row_size = %d, Jsize = %d'));
    %    cn(rst(i):rst(i+1)-1) = J;
    %    val(rst(i):rst(i+1)-1) = AA;

    A_ell.ja(i, 1:rl(i)) = J(:);
    A_ell.aa(i, 1:rl(i)) = AA(:);
end