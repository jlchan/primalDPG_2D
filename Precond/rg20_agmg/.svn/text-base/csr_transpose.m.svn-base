function At = csr_transpose(A)
% function At = csr_transpose(A)
% purpose : returns transpose of a give csr matrix


% convert to matlab format
A_mat = csr2mat(A);

% transpose
At_mat = A_mat';

% convert back to csr format
At = mat2csr(At_mat);
