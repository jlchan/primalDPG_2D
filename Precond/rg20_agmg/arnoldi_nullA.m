 function  [V,H, kact]=arnoldi_nullA(A, nullA, v0, k, reorth)
%
%-----------------------------------------------------------
%  k-step Arnoldi:    A V(:,1:k) = V H, V(:,1)=v0/norm(v0).
%-----------------------------------------------------------
%
%  Input:
%
%      A      (n-by-n) the matrix
%      v0     n-vector
%      k      # of Arnoldi steps requested
%    reorth   (optional) set to 1 for reorthogonalization, (default)
%             set to any other value to switch it off
%
%  Output:
%
%      kact   actual Arnoldi steps taken
%
%             -----  kact=k -------
%      V      n-by-(k+1)  Arnoldi vectors
%      H      (k+1)-by-k
%             -----  kact=j<k -------
%      V      n-by-j  Arnoldi vectors
%      H      j-by-j
%
% (c) Ren-Cang Li, rcli@uta.edu,  06/16/07

if nargin < 3,
    disp('not enough inputs');
    return
elseif nargin == 3
    reorth=1;
end
[m,n]=size(A);
if m ~= n,
    disp('the sizes of input matrix incorrect');
    return
end
%
V = zeros(n,k+1); H = zeros(k+1,k);

nullA = nullA/norm(nullA);

Dinv = diag(A) + nullA.*nullA;
Dinv = 1./Dinv;

nrm2 = norm(v0);
if nrm2 == 0.0,
   disp('arnoldi: input v0 is a zero vector');
   return
end
tol=n*eps;

V(:,1)=v0/nrm2;
for j=1:k,
    vh = Dinv.*(A*V(:,j) + nullA*(nullA'*V(:,j)));
    nrmvh=norm(vh);
%   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%   by MGS
    for i=1:j,
        hij=( V(:,i) )' * vh;
        vh = vh - hij*V(:,i);
        H(i,j) = hij;
    end
    if reorth == 1,
       for i=1:j,
           tmp=( V(:,i) )' * vh;
           vh = vh - tmp*V(:,i);
           H(i,j)=H(i,j)+tmp;
       end
    end
%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    H(j+1,j)=norm(vh);
    V(:,j+1)=vh/H(j+1,j);
    if H(j+1,j) <= tol*nrmvh,
       disp(strcat('termination at step:', num2str(j)));
       kact = j; return;
    end
 end
 kact = k;