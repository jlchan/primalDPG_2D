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
Ntrial = [1:8];

% grids = {'Maxwell0125.neu'};
% Ntrial = 4;
% plotFlag = 1;

Ntest = Ntrial+2;
Nflux = Ntrial;
resVecs = {};
for i = 1:length(grids)
    for j = 1:length(Ntrial)
        [A, b, nU, nM, Np, Rr bK] = primalDPG_poisson(grids{i},Ntrial(j),Ntest(j),Nflux(j),plotFlag);
        if ~plotFlag
            Av = A(1:nU,1:nU); Af = A(nU+(1:nM),nU+(1:nM));Avf = A(1:nU,nU+(1:nM));
            % Pre = cholinc(A,1e-2); %Pre'*Pre
            % Pre = chol(A);
            % Pfun = @(x) Pre\(( Pre')\x);            
            keyboard
            Blk = [Av zeros(size(Avf));zeros(size(Avf')) Af]; 
%             Pre = @(x) Rr*(Kbd\(Rr'*x)); % edge overlap only
            Pre = @(x) Blk\x; % "ideal" block diag precond
            
            % [U, flag, relres, iter, resvec] = pcg(A,b,1e-7,250,Pre',Pre);
            [U, flag, relres, iter, resvec] = pcg(A,b,1e-8,250,@(x) Pre(x));
            resVecs{i,j} = resvec;
        end
    end
end

[m n] = size(resVecs);
C = hsv(n*m);
for i = 1:m
    for j = 1:n
        semilogy(resVecs{i,j},'.-','color',C((i-1)*n+j,:));
        hold on
    end
end