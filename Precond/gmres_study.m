function gmres_study

% Initializes paths
% addpath .
% addpath ../Codes1D
% addpath ../Codes2D
% addpath ../CodesDPG
% addpath ../ServiceRoutines
% addpath ../Grid/
% addpath ../Grid/Maxwell2D
% addpath ../Grid/Other
addpath('rg20_agmg')

plotFlag = 0;
grids = {'Maxwell05.neu','Maxwell025.neu','Maxwell0125.neu'};
grids = {'Maxwell05.neu','Maxwell025.neu'};
% grids = {'Maxwell00625.neu'};
% grids = {'Maxwell025.neu'};
Ntrial = [2 6];
NpTrials = (Ntrial+1).*(Ntrial+2)/2;

% grids = {'Maxwell1.neu'};
% Ntrial = 4;
% plotFlag = 1;

Ntest = Ntrial+4;
Nflux = Ntrial;
resVecs = {};
legendVec = {};
for i = 1:length(grids)
    for j = 1:length(Ntrial)
        NpTrial = NpTrials(j);
        % for poisson, set b = 0, eps = 1
        b = 0; eps = 1;
%         b = 1; eps = 0e-6;
        [A, b, nU, nM, Np, Rp, Irp, M, fpairs] = primalDPG_confusion(grids{i},Ntrial(j),Ntest(j),Nflux(j),plotFlag,b,eps);
        if plotFlag            
            keyboard
            S = Af-Avf'*(Av\Avf);
        end
        ui = 1:nU; mi = nU + (1:nM);
        Av = A(ui,ui); Af = A(mi,mi);Avf = A(ui,mi);
        %             
        %             Sx = @(x) Af*x-Avf'*(Av\(Avf*x));
        precondFlag = 'oas';
        switch precondFlag
            case 'ideal'
                % build true preconditioner
                AvPre = @(x) Av\x;
                SPre = @(x) Af\x; % ignore Schur complement...
            case 'oas'                
                AvPre = buildOAS_CG(Rp,Av,Ntrial(j));%
%                 AvPre = @(x) Av\x;
%                 SPre = @(x) Af\x; % ignore Schur complement...
                SPre = buildOAS_mortar(Af,Nflux(j),fpairs);
%                 R = ichol(Af);
%                 SPre = @(x)R\(R'\x);
%                 levelsAf = agmg_setup(Af);
%                 SPre = @(x) agmg_solve(levelsAf,x,5,1e-3);

%                 AvPre = @(b) fgmres(Av,b,1e-16, 2, buildOAS_CG(Rp,Av,Ntrial(j)));
%                 SPre  = @(b) fpcg(Af,b,1e-6, 2, buildOAS_mortar(Af,Nflux(j)));
%                 bu = rand(size(Av,1),1); bf = rand(size(Af,1),1);
%                 [u, flag, relres, iter, resvecu] = fpcg(Av,bu,1e-7,100,@(x) AvPre(x));
%                 [f, flag, relres, iter, resvecf] = fpcg(Af,bf,1e-7,100,@(x) SPre(x));
%                 keyboard
            case 'agmg'
                % use AMG preconditioner
                levelsAv = agmg_setup(Av);
                %                     [x flag relres iter resvec] = agmg_solve(levels, b, maxiter, tol);
                AvPre = @(x) agmg_solve(levelsAv,x,50,1e-2);
                %                     AvPre = @(x) Av\x;
                levelsAf = agmg_setup(Af);
                SPre = @(x) agmg_solve(levelsAf,x,50,1e-3);
                %                     SPre = @(x) S\x;
                %                     SPre = @(x) Af\x;                               
        end
        
        %% Define fixed point iteration precond
        iters = 1;
        Pre = @(x) fixedPoint(x,AvPre,Avf,SPre,ui,mi,iters);
        
        %% GMRES        
        [U, flag, relres, iter, resvec] = fgmres(A,b,1e-6,100,@(x) Pre(x));
        resVecs{i,j} = resvec;
    end
end

[m n] = size(resVecs);
C = hsv(n*m);
% figure
gcf;semilogy(0,0)
for i = 1:m
    for j = 1:n
        semilogy(resVecs{i,j},'.-','color',C((i-1)*n+j,:));
        hold on
        legendVec{(i-1)*n+j} = ['grid ' num2str(i) ', N = ', num2str(Ntrial(j))];
    end
end
legend(legendVec);
xlabel('Iteration count')
ylabel('Residual')
title([precondFlag ' preconditioner'])

function x = fixedPoint(b,Ainv,B,Cinv,I1,I2,iters)

f = b(I1);
g = b(I2);
% u1 = zeros(length(I1),1);
% u2 = zeros(length(I2),1);
u1 = f;
u2 = g;
for i = 1:iters
    tmp = (b - [B*u2;B'*u1]);
    x = [Ainv(tmp(I1)); Cinv(tmp(I2))];
    u1 = x(I1);u2 = x(I2);
end
x(I1) = u1;
x(I2) = u2;