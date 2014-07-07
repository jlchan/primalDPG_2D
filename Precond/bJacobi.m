% apply block jacobi iteration
function [U r] = bJacobi(A,b,nU,nM,nIter)

if nargout>1
    Uex = A\b;
    uex = Uex(1:nU);
end

uh = zeros(nM,1);

I1 = 1:nU; I2 = nU + (1:nM);

A1 = A(I1,I1); B = A(I1,I2); C = A(I2,I2);
b1 = b(I1); b2 = b(I2);
% levelsC = agmg_setup(C);
r = zeros(nIter,1);
for k = 1:nIter
    %         Mr = Rr*M*Rr';
    %         D = blkdiag(A1,C);
    %         u =  (Rr*M*Rr')\(b1-B*uh);
    u =  A1\(b1-B*uh);
    uh = C\(b2-B'*u);
%     keyboard
%     [uh flag relres iter resvec] = agmg_solve(levelsC, b2-B'*u, 50, 1e-7);
    if (nargout>1)
        r(k) = norm(u-uex);
    end
    
    %         u = Rr'*u;
    %         Nplot = Ntrial; [xu,yu] = Nodes2D(Nplot);
    %         Nplot = 25; [xu,yu] = EquiNodes2D(Nplot);
    %         [ru, su] = xytors(xu,yu);
    %         Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
    %         xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
    %         yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
    %         color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.');
    %         title(['k = ' num2str(k)])
    %         pause
end
U = zeros(size(A,1),1);
U(I1) = u;
U(I2) = uh;

