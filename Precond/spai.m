% simple (and probably slow) sparse approx inverse preconditioner 
% from Grote/Huckle

function iA = spai(A,tol)

[m, n] = size(A);
iA = sparse(m,n);
for k = 1:n % assume all diags nonzero   
    J = k; % start initial sparsity = diagonal
    e = zeros(m,1); e(k) = 1;
    mk = A(:,k)\e;    
    r = A(:,k)*mk - e;    
    r2 = r'*r;
    % add more indices 
    while r2 > tol*tol        
        [~, Jh] = find(A(abs(r)> 10*eps,:));                
        Jh = setdiff(unique(Jh),J); 
        
        % find new indices by optimizing over residual reduction
        rTAej = r'*A(:,Jh);  
        rho = zeros(length(Jh),1);
        for j = 1:length(Jh)
            rho(j) = r2 - (rTAej(j)/norm(A(:,Jh(j))))^2;
        end                        
        newJ = Jh(abs(rho-min(rho)) < 1e-8);        
        J = [J;newJ(:)];
        J = unique(J);
        
        mk = A(:,J)\e;
        r = A(:,J)*mk - e; 
        r2 = r'*r;    
    end    
    iA(J,k) = mk;
end

