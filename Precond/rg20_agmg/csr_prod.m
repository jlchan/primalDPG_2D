function AB = csr_prod(A, B)
% function AB = csr_prod(A, B)
% purpose : returns product of two csr matrices A and B

% convert A and B to matlab format
A_mat = csr2mat(A);
B_mat = csr2mat(B);

% do the product with matlab format
AB_mat = A_mat*B_mat;

% convert back to csr format
AB = mat2csr(AB_mat);
