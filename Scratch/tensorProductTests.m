% experiments w/tensor products
function tensorProductTests
Globals2D
[Nv, VX, VY, K, EToV] = QuadMesh2D(1);
kL = zeros(10,1);
kxy = kL;
for N = 1:10
    % N = 10; % when N = even, Nf = N-1, fails?
    
    % Initialize solver and construct grid and metric
    StartUp2D
    
    Dr = full(Dr);Ds = full(Ds);
    M = MassMatrix;
    Np1D = N+1;
    r1D = x(1:Np1D);
    V1D = Vandermonde1D(N,r1D);
    M1D = inv(V1D*V1D');
    D1D = Dmatrix1D(N,r1D,V1D);
    
    K1D = D1D'*M1D*D1D;    
    Kx = kron(M1D,K1D); % x is stored in the second half
    Ky = kron(K1D,M1D); % x is stored in the second half
        
    Kh = Dr'*MassMatrix*Dr + Ds'*MassMatrix*Ds;
    f = ones(Np1D^2,1);
    b = M*f;
    
    L = D1D^2;
    % L(1,:) = 0;L(:,1) = 0;L(1,1) = 1;
    % L(Np1D,:) = 0;L(:,Np1D) = 0;L(Np1D,Np1D) = 1;
    Lap = kron(L,eye(Np1D)) + kron(eye(Np1D),L); 
    
    I = reshape(1:Np1D^2,Np1D,Np1D);
    i = I(1,:);
    Lap = applyZeroBCs(Lap,i);f(i) = 0;
    Kh = applyZeroBCs(Kh,i);b(i) = 0;
    Kx(i,:) = 0; Kx(:,i) = 0;Kx(i,i) = .5*eye(length(i));
    Ky(i,:) = 0; Ky(:,i) = 0;Ky(i,i) = .5*eye(length(i));
    i = I(Np1D,:);
    Lap = applyZeroBCs(Lap,i);f(i) = 0;
    Kh = applyZeroBCs(Kh,i);b(i) = 0;
    Kx(i,:) = 0; Kx(:,i) = 0;Kx(i,i) = .5*eye(length(i));
    Ky(i,:) = 0; Ky(:,i) = 0;Ky(i,i) = .5*eye(length(i));
    i = I(:,1);
    Lap = applyZeroBCs(Lap,i);f(i) = 0;
    Kh = applyZeroBCs(Kh,i);b(i) = 0;
    Kx(i,:) = 0; Kx(:,i) = 0;Kx(i,i) = .5*eye(length(i));
    Ky(i,:) = 0; Ky(:,i) = 0;Ky(i,i) = .5*eye(length(i));
    i = I(:,Np1D);
    Lap = applyZeroBCs(Lap,i);f(i) = 0;
    Kh = applyZeroBCs(Kh,i);b(i) = 0;
    Kx(i,:) = 0; Kx(:,i) = 0;Kx(i,i) = .5*eye(length(i));
    Ky(i,:) = 0; Ky(:,i) = 0;Ky(i,i) = .5*eye(length(i));
%     plotSol(-Lap\f,50);
    % plotSol(Kh\b,50);
    % plotSol(Kx\b + Ky\b,50);
    
    kL(N) = cond(-Lap\Kh);
    kxy(N) = cond(Kx\Kh + Ky\Kh);
end

semilogy(kL,'r-');hold on
semilogy(kxy,'b-')

function A = applyZeroBCs(A,i)

A(i,:) = 0; A(:,i) = 0;A(i,i) = eye(length(i));
