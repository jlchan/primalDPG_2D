% gets 1D DPG ultra-weak operator for poisson on a uniform grid of K
% elements, with test order NT and trial order N
% function UW_DPG_1D

Kmax = 8;
Nvec = [4];
Kvec = [2.^(4:Kmax)];

c = cell(length(Kvec),1);
iter = c;
for K = Kvec    
    c{K} = zeros(4,Nmax);
    iter{K} = zeros(4,Nmax);
    for N = Nvec
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
%         c{K}(1,N) = condest(A);
%         c{K}(2,N) = condest(P1*A);
%         c{K}(3,N) = condest(P2*A);
%         c{K}(4,N) = condest(P0(A) + P1*A);
%         c{K}(5,N) = condest(P0(A) + P2*A);
        
%         [~,~,~,it1,resvec1] = pcg(A,b,1e-10,100);
%         [~,~,~,it2,resvec2] = pcg(A,b,1e-10,100,@(x) P1*x);
%         [~,~,~,it3,resvec3] = pcg(A,b,1e-10,100,@(x) P2*x);
        it1 = nan; it2 = nan; it3 = nan;
        [~,~,~,it4,resvec4] = pcg(A,b,1e-10,100,@(x) P0(x) + P1*x);
        [~,~,~,it5,resvec5] = pcg(A,b,1e-10,100,@(x) P0(x) + P2*x);
        iter{K}(1,N) = it1;
        iter{K}(2,N) = it2;
        iter{K}(3,N) = it3;
        iter{K}(4,N) = it4;
        iter{K}(5,N) = it5;
    end    
    semilogy(resvec4,'.-');hold on;
    semilogy(resvec5,'r.-');hold on;
end
break

       
%%
% for K = 2.^(2:Kmax)
%     figure
% semilogy(c{K}','.-');
%     legend('No pre','0 overlap','1 overlap','1 overlap + P0 coarse grid')
%     title(sprintf('Num elements = %i',K))
%     xlabel('polynomial order')
%     ylabel('Condition number')
% end        
        
for N = Nvec
    c1 =[]; c2 =[]; c3 =[]; c4 =[];  c5 = [];
    it1 = []; it2 = []; it3 = []; it4 = []; it5 = [];
    for K = 2.^(2:Kmax)    
%         c1 =[c1 c{K}(1,N)];
%         c2 =[c2 c{K}(2,N)];
%         c3 =[c3 c{K}(3,N)];
%         c4 =[c4 c{K}(4,N)];
%         c5 =[c5 c{K}(5,N)];
        
        it1 = [it1 iter{K}(1,N)];
        it2 = [it2 iter{K}(2,N)];
        it3 = [it3 iter{K}(3,N)];
        it4 = [it4 iter{K}(4,N)];        
        it5 = [it5 iter{K}(5,N)];        
    end
%     figure    
%     loglog(2.^(2:Kmax),c1,'r.-');hold on
%     semilogy(2.^(2:Kmax),c2,'g.-')
%     semilogy(2.^(2:Kmax),c3,'b.-')
%     semilogy(2.^(2:Kmax),c4,'k.-')
%     semilogy(2.^(2:Kmax),c5,'m.-')
%     xlabel('Number of elements')
%     ylabel('Condition number')
%     legend('No pre','0 overlap','1 overlap','0 overlap + P0 coarse grid','1 overlap + P0 coarse grid')
%     title(sprintf('Polynomial order = %i',N))
    
    figure
    loglog(2.^(2:Kmax),it1,'r.-');hold on
    loglog(2.^(2:Kmax),it2,'g.-');
    loglog(2.^(2:Kmax),it3,'b.-');
    loglog(2.^(2:Kmax),it4,'k.-');
    loglog(2.^(2:Kmax),it5,'m.-');
    xlabel('Number of elements')
    ylabel('Iterations to convergence')

    legend('No pre','0 overlap','1 overlap','0 overlap + P0 coarse grid','1 overlap + P0 coarse grid')
    title(sprintf('Polynomial order = %i',N))
end        
