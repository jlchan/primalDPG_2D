clear

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
grids = {'Maxwell1.neu','Maxwell05.neu','Maxwell025.neu','Maxwell0125.neu'};
grids = {'Maxwell1.neu','Maxwell05.neu','Maxwell025.neu'};
Ntrial = [1:2];

grids = {'Maxwell1.neu'};
Ntrial = 4;
plotFlag = 1;

Ntest = Ntrial+2;
Nflux = Ntrial;
resVecs = {};
legendVec = {};
for i = 1:length(grids)
    for j = 1:length(Ntrial)
        [A, b, nU, nM, Np, bK, Rp, Irp, M] = primalDPG_poisson(grids{i},Ntrial(j),Ntest(j),Nflux(j),plotFlag);
        if ~plotFlag
            Av = A(1:nU,1:nU); Af = A(nU+(1:nM),nU+(1:nM));Avf = A(1:nU,nU+(1:nM));            
            
            R = [speye(nU) Av\Avf; zeros(nM,nU) speye(nM)];
            Rinv = [speye(nU) -Av\Avf; zeros(nM,nU) speye(nM)];            
            P = [Av zeros(nU,nM);zeros(nM,nU) Af - Avf'*(Av\Avf)];            
            
            
%             S = Af;
%             Blk = [Av zeros(size(Avf));zeros(size(Avf')) S];
%             Pre = @(x) Blk\x;
            
            ui = 1:nU; mi = nU + (1:nM);            
            AvPre = @(x) Av\x;
%             AvPre = @(x) Rp*(bK\(Rp'*x));
            RTinv = @(x) [x(ui); x(mi) - Avf'*AvPre(x(ui))];
            Rinv = @(x) [x(ui)-AvPre(Avf*x(mi)); x(mi)];
                        
            SPre = @(x) Af\x; %Af*x-Avf'*AvPre(Avf*x);
            P = @(x) [AvPre(x(1:nU)); SPre(x(mi))];          %[Av zeros(nU,nM);zeros(nM,nU) Af - Avf'*AvPre(Avf)];
            Pre = @(x) Rinv(P(RTinv(x)));
%             Pre = @(x) P(x);
%             keyboard
                        
            % [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,250,Pre',Pre);
            [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,50,@(x) Pre(x));
            resVecs{i,j} = resvec;            
        end
    end
end

[m n] = size(resVecs);
C = hsv(n*m);
figure
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
% title('Block additive preconditioner')
% title('Block ideal preconditioner')