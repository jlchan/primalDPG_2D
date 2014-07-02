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
% grids = {'Maxwell05.neu','Maxwell025.neu','Maxwell0125.neu'};
grids = {'Maxwell05.neu','Maxwell025.neu'};
% grids = {'Maxwell00625.neu'};
grids = {'Maxwell025.neu'};
Ntrial = [2:6];
NpTrials = (Ntrial+1).*(Ntrial+2)/2;

% grids = {'Maxwell1.neu'};
% Ntrial = 4;
% plotFlag = 1;

Ntest = Ntrial+2;
Nflux = Ntrial;
resVecs = {};
legendVec = {};
for i = 1:length(grids)
    for j = 1:length(Ntrial)       
        NpTrial = NpTrials(j);
%                 [A, b, nU, nM, Np, Rp, Irp, M] = primalDPG_poisson(grids{i},Ntrial(j),Ntest(j),Nflux(j),plotFlag);
        [A, b, nU, nM, Np, Rp, Irp, M] = primalDPG_confusion(grids{i},Ntrial(j),Ntest(j),Nflux(j),plotFlag); 
        if ~plotFlag
            
            Av = A(1:nU,1:nU); Af = A(nU+(1:nM),nU+(1:nM));Avf = A(1:nU,nU+(1:nM));            
%             S = Af-Avf'*(Av\Avf);
            Sx = @(x) Af*x-Avf'*(Av\(Avf*x));
            precondFlag = 'agmg';
            switch precondFlag
                case 'ideal'
                    % build true preconditioner
                    AvPre = @(x) Av\x;
                    SPre = @(x) Af\x; % ignore Schur complement...
%                     SPre = @(x) S\x;
%                     SPre = @(x) Mfaces\x;
                case 'oas'
                    % build OAS preconditioner
                    RAv = Rp'*Av*Rp; % expand out to local dofs
                    K = max(size(Rp))/NpTrial;
                    OAS = {};
                    for k = 1:K
                        inds = NpTrial*(k-1) + (1:NpTrial);
                        OAS{k} = RAv(inds,inds);
                    end
                    OAS = blkdiag(OAS{:});
                    AvPre = @(x) Rp*(OAS\(Rp'*x));
                    SPre = @(x) S\x; 
                case 'agmg'
                    % use AMG preconditioner
                    levelsAv = agmg_setup(Av);
                    levelsAf = agmg_setup(Af);
                    tolAv = 5e-3;
                    maxiterAv = 50;
                    tolAf = 1e-3;
                    maxiterAf = 50;
                    %             [x flag relres iter resvec] = agmg_solve(levels, b, maxiter, tol);
%                     AvPre = @(x) agmg_solve(levelsAv,x,maxiterAv,tolAv);
                    AvPre = @(x) Av\x;
%                     SPre = @(x) agmg_solve(levelsAf,x,maxiterAf,tolAf);                    
%                     SPre = @(x) S\x;
                    SPre = @(x) Af\x;
                case 'fem-agmg'
                    
            end
            
            % build block cholesky bits
            ui = 1:nU; mi = nU + (1:nM);
            RTinv = @(x) [x(ui); x(mi) - Avf'*AvPre(x(ui))];
            Rinv = @(x) [x(ui)-AvPre(Avf*x(mi)); x(mi)];
            
            %
            PBlockDiag = @(x) [AvPre(x(1:nU)); SPre(x(mi))];          %[Av zeros(nU,nM);zeros(nM,nU) Af - Avf'*AvPre(Avf)];
            Pre = @(x) Rinv(PBlockDiag(RTinv(x)));
%             Pre = @(x) PBlockDiag(x);
            %             keyboard
            
            %% PCGGGGGGGGGGGGGG
            % [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,250,Pre',Pre);
            [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,50,@(x) Pre(x));
            resVecs{i,j} = resvec;            
        end
    end
end

[m n] = size(resVecs);
C = hsv(n*m);
figure
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
