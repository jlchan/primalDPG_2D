% gets 1D DPG ultra-weak operator for poisson on a uniform grid of K
% elements, with test order NT and trial order N
% function UW_DPG_1D
clear
for N = [2 4 8]
    disp(sprintf('Poly order %i',N))
    for K = [2 4 8 16 32 64]    
        % N = 1;        
        [A b ids] = UW_DPG_1D(N,K);
        
        % preconditioners
        P1 = zeros(size(A,1),size(A,2));
        P2 = zeros(size(A,1),size(A,2));
        for e = 1:K
            % zero overlap
            [U S V] = svd(full(A(ids{e},ids{e})));
            P1(ids{e},ids{e}) = P1(ids{e},ids{e}) + V*diag(1./diag(S))*U';
            
            % one element overlap
            if e==1
                idK = unique([ids{e} ids{e+1}]);
            elseif e==K
                idK = unique([ids{e-1} ids{e}]);
            else
                idK = unique([ids{e-1} ids{e} ids{e+1}]);
            end            
            [U S V] = svd(full(A(idK,idK)));
            P2(idK,idK) = P2(idK,idK) + V*diag(1./diag(S))*U';
        end
        
        % % P1 coarse grid instead
%         r0 = JacobiGL(0,0,1); I0 = Vandermonde1D(1,r) * inv(Vandermonde1D(1,r0));        
        r = JacobiGL(0,0,N);
        I0 = ones(size(r));
        
        I0 = kron(speye(K),I0);
        I0 = blkdiag(I0,I0,eye(2*(K+1)));
        P0 = @(x) I0*((I0'*A*I0)\(I0'*x));
        c1=cond(A);
        c2=cond(P1*A);
        c3=cond(P2*A);
        c4=cond(P0(A) + P2*A);
        c5=cond(P0(A) + P1*A);
        [~,~,~,it1] = pcg(A,b,1e-10,250);
        [~,~,~,it2] = pcg(A,b,1e-10,250,@(x) P1*x);
        [~,~,~,it3] = pcg(A,b,1e-10,250,@(x) P2*x);
        [~,~,~,it4] = pcg(A,b,1e-10,250,@(x) P0(x) + P2*x);       
        [~,~,~,it5] = pcg(A,b,1e-10,250,@(x) P0(x) + P1*x);
        %disp(sprintf('%i cells: condnums = %1.1d, %1.1d, %1.1d, and iterations = %i, %i, %i, %i',K,c2,c3,c4,it2,it3,it4))        
        disp(sprintf('%i cells: condnums = %1.1d, %1.1d, and iterations = %i, %i',K,c4,c5,it4,it5))
    end    
    disp(sprintf('\n'))
end
       

