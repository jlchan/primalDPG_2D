function A_csr = mat2csr(A)
% function A_csr = mat2csr(A)
% purpose : converts a matrix stored in matlab format to csr format
n = size(A,1);
m = size(A,2);
nz = nnz(A);

A_csr = csr(n,m,nz);

ia = A_csr.ia;
ja = A_csr.ja;
aa = A_csr.aa;

row_counters  = zeros(n,1);
counter = 1;

ia(1) = counter;
for i=1:n
    [I J AA] = find(A(i,:));

    I = 0*I+i;
    % length(I)
    % length(J)
    % length(AA)

    % make sure that diagonal is first
    if(n == m)
        temp = abs(I(:) - J(:));
        [vals ids] = sort(temp);

        if(length(I) ~= length(ids))
            display('something is wrong');
            return;
        end

        I = I(ids);
        J = J(ids);
        AA = AA(ids);

        % if(J(1) ~= i)
        %     i
        %     display('error in converting to csr...');
        %     return;
        % end

    end

    row_size = length(I);
    if(row_size > 0)
        ja(counter:counter+row_size-1) = J(:);
        aa(counter:counter+row_size-1) = AA(:);
        counter = counter + length(I);
    end
    ia(i+1) = counter;
end

A_csr.ia = ia;
A_csr.ja = ja;
A_csr.aa = aa;

return