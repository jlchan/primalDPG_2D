function Mhandle = buildOAS_CG(R,A,Norder)

Globals2D;

if nargin<3
    Norder = N;
end

% build OAS preconditioner algebraically:
Norderp = (Norder+1)*(Norder+2)/2;
[r c] = find(R); r = reshape(r,Norderp,K);
Ak = cell(K,1);Aki = cell(K,1);
for k = 1:K
%     nbr = unique([k EToE(k,:)]); % 1 elem overlap
    nbr = k; % no overlap
    inds = unique(r(:,nbr));
    Aki{k} = inds;
    Ak{k} = eye(numel(inds));%A(inds,inds);
end

% coarse grid solver
[Rc1 Ir1 vmapBT1 xr1 yr1] = pRestrictCG(Norder,1); % interp down
Rp1 = Rc1*Ir1';
Rs = spdiag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
R1 = Rs*Rp1'; % interp down to P1
P1 = @(b) R1*((R1'*A*R1)\(R1'*b)); % should be able to "ignore" BC imposition - pcg acts on residual e(vmapBT) = 0

Mhandle = @(x) P1(x) + OAS(x,Ak,Aki); % coarse solver + OAS solver

