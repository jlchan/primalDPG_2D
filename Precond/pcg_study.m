% Initializes paths
addpath .
addpath ../Codes1D
addpath ../Codes2D
addpath ../CodesDPG
addpath ../ServiceRoutines
addpath ../Grid/
addpath ../Grid/Maxwell2D
addpath ../Grid/Other

plotFlag = 0;
grids = {'Maxwell05.neu','Maxwell025.neu','Maxwell0125.neu'};
grids = {'Maxwell05.neu','Maxwell025.neu'};
% grids = {'Maxwell00625.neu'};
% grids = {'Maxwell05.neu'};
Ntrial = [2 4];
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
                %                     SPre = @(x) S\x;
                %                     SPre = @(x) Mfaces\x;
            case 'oas'
                %                     % build OAS preconditioner
                %                     RAv = Rp'*Av*Rp; % expand out to local dofs
                %                     K = max(size(Rp))/NpTrial;
                %                     OAS = {};
                %                     for k = 1:K
                %                         inds = NpTrial*(k-1) + (1:NpTrial);
                %                         OAS{k} = RAv(inds,inds);
                %                     end
                %                     OAS = blkdiag(OAS{:});
                %                     AvPre = @(x) Rp*(OAS\(Rp'*x));
                %                     SPre = @(x) S\x;
                
                AvPre = buildOAS_CG(Rp,Av,Ntrial(j));%
%                 AvPre = @(x) Av\x;
%                                     SPre = @(x) Af\x; % ignore Schur complement...
                SPre = buildOAS_mortar(Af,Nflux(j),fpairs);
% S = Af-Avf'*(Av\Avf);SPre = buildOAS_mortar(S,Nflux(j),fpairs);                
%                 AvPre = @(b) fpcg(Av,b,1e-16, 2, buildOAS_CG(Rp,Av,Ntrial(j)));
%                 SPre  = @(b) fpcg(Af,b,1e-6, 2, buildOAS_mortar(Af,Nflux(j)));
                bu = rand(size(Av,1),1); bf = rand(size(Af,1),1);
                [u, flag, relres, iter, resvecu] = fpcg(Av,bu,1e-7,100,@(x) AvPre(x));
                [f, flag, relres, iter, resvecf] = fpcg(Af,bf,1e-7,100,@(x) SPre(x));
                keyboard
            case 'agmg'
                % use AMG preconditioner
                levelsAv = agmg_setup(Av);
                levelsAf = agmg_setup(Af);
                tolAv = 5e-3;
                maxiterAv = 50;
                tolAf = 1e-3;
                maxiterAf = 50;
                %                     [x flag relres iter resvec] = agmg_solve(levels, b, maxiter, tol);
                AvPre = @(x) agmg_solve(levelsAv,x,maxiterAv,tolAv);
                %                     AvPre = @(x) Av\x;
                SPre = @(x) agmg_solve(levelsAf,x,maxiterAf,tolAf);
                %                     SPre = @(x) S\x;
                %                     SPre = @(x) Af\x;
                
            case 'fem-agmg'
                
        end
        
        % build block cholesky bits
        RTinv = @(x) [x(ui); x(mi) - Avf'*AvPre(x(ui))];
        Rinv = @(x) [x(ui)-AvPre(Avf*x(mi)); x(mi)];
        
        PBlockDiag = @(x) [AvPre(x(1:nU)); SPre(x(mi))];          %[Av zeros(nU,nM);zeros(nM,nU) Af - Avf'*AvPre(Avf)];
        Pre = @(x) Rinv(PBlockDiag(RTinv(x)));
        %             Pre = @(x) PBlockDiag(x);
        %             keyboard
        %             nIter = 3;
        %             Pre = @(x) bJacobi(A,x,nU,nM,nIter);
        
        %% PCGGGGGGGGGGGGGG
        % [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,250,Pre',Pre);
        [U, flag, relres, iter, resvec] = fgmres(A,b,1e-6,100,@(x) Pre(x));
        resVecs{i,j} = resvec;
    end
end

[m n] = size(resVecs);
C = hsv(n*m);
% figure
semilogy(0,0)
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
