clear
grids = {'Maxwell05.neu','Maxwell025.neu','Maxwell0125.neu'};
addpath ../CodesDPG
addpath rg20_agmg
mesh = grids{2};
% mesh = 'Maxwell1.neu';

Globals2D;
N = 6;
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);
StartUp2D;

% build op
[M, Dx, Dy] = getBlockOps();
a = 0;
AK = a*M + Dx'*M*Dx + Dy'*M*Dy;
[R, vmapBT] = getCGRestriction();
A = R*AK*R';
f = ones(Np*K,1);
b = R*M*f;
    
u0 = zeros(size(A,1),1);

% Dirichlet BC on the left
left = x(vmapB) < -1+NODETOL;
u0(vmapBT) = 0;
vmapBT(~left) = []; % remove left BCs for Neumann
b = b - A*u0;
b(vmapBT) = u0(vmapBT);
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));
    
u = R'*(A\b);

Nplot = 25;
[xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.')

preFlag = '';
preFlag = 'nodalfem';
% preFlag = 'oas';
switch preFlag
    case 'oas'
        % OAS preconditioner + P1 coarse
        [Rp, Irp, vmapBTr, xr, yr] = pRestrictCG(1); % restrict test to trial space
        Rr = Rp*Irp'; % interp down to linears
        A1 = Rr*AK*Rr';

        left = xr < -1+NODETOL;
        vmapBTr(~left) = []; % remove left BCs for Neumann        
        A1(vmapBTr,:) = 0; A1(:,vmapBTr) = 0;
        A1(vmapBTr,vmapBTr) = speye(length(vmapBTr));
        
        % build OAS preconditioner
        RA = R'*A*R; % expand out assembled mat to local dofs
        OAS = {};
        for k = 1:K
            inds = Np*(k-1) + (1:Np);
            OAS{k} = RA(inds,inds);
        end
        OAS = blkdiag(OAS{:});
                
        % scale by # of connections 
        Pre = @(x) diag(1./sum(R,2))*R*(OAS\(R'*x)) + A1;
        keyboard
    case 'nodalfem'
        [r c] = find(R); [ru i] = unique(r);
        xu = x(c(i)); yu = y(c(i));   % get unique nodes
        tri = delaunay(xu,yu);        
        
        % make new mesh from delaunay
        EToV = tri; 
        K = size(tri,1);  Nv = length(xu);
        VX = xu(:)';VY = yu(:)';
        oldN = N;N = 1; StartUp2D;
        [M, Dx, Dy] = getBlockOps();
        AK1 = a*M + Dx'*M*Dx + Dy'*M*Dy;
        [R1, vmapBT] = getCGRestriction();
        A1 = R1*AK1*R1';
        left = x(vmapB) < -1+NODETOL;
        vmapBT(~left) = []; % remove left BCs for Neumann
        b1 = R1*M*ones(Np*K,1);
        b1(vmapBT) = 0;
        A1(vmapBT,:) = 0; A1(:,vmapBT) = 0;
        A1(vmapBT,vmapBT) = speye(length(vmapBT));
        %         u = R'*(A1\b);
        
        %         plot(A\b);hold on;plot(A1\b,'r')
        levels = agmg_setup(A);
        levels1 = agmg_setup(A1);
        maxiter = 5;
        tol = 1e-4;
        
        [x flag relres iter resvec1] = agmg_solve(levels1, b, maxiter, tol);
        [UD flag relres iter resvec] = agmg_solve(levels, b, maxiter, tol);
        
        semilogy(resvec,'.-');hold on;
        semilogy(resvec1,'r.-')
        %         legend('Inner iteration cost of AGMG')
        legend(['N = ' num2str(oldN)],'P1 FEM')
                
        Pre1 = @(x) agmg_solve(levels1,x,maxiter,tol);
        Pre = @(x) agmg_solve(levels,x,maxiter,tol);
end
levels = agmg_setup(A);
[UN, flag, relres, iter, resvecN] = agmg_solve(levels,b,100,1e-4); % using gmres
% [UNP, flag, relres, iter, resvecNp] = pcg(A,b,1e-4,100,@(x) Pre(x));flag
[U1, flag, relres, iter, resvec1] = pcg(A,b,1e-4,100,@(x) Pre1(x));flag

figure
semilogy(resvecN,'.-');hold on;
% semilogy(resvecNp,'k.-'); 
semilogy(resvec1,'r.-')
legend(['Direct AGMG, N = ' num2str(oldN)],'AGMG inner pre','P1 FEM AGMG preconditioner')
break
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh);N = oldN; StartUp2D; % reset to original mesh
% u = R'*UN;  
% UD = A\b;
u = R'*UN;

Nplot = 25;
[xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.')
