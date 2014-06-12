% returns R such that R*Ak*R' gives back a CG matrix from the
% block diagonal matrix Ak.
% other return values: vmapBT, the boundary indices in the reduced system.  
function [R vmapBT] = getCGRestriction()

Globals2D

% unique vertex nodes: search face nodes for vertex matches
Npv = numel(VX);vmapT = zeros(Np*K,1);
for i = 1:Npv
    inds = find(abs(x(vmapM)-VX(i))<1e-14 & abs(y(vmapM)-VY(i))<1e-14);
    vmapT(vmapM(inds)) = i;
end

% unique face dofs - find non-vertex face dofs
nonVerts = ~ismember(vmapM,find(vmapT)); 
vmapMF = vmapM(nonVerts); vmapPF = vmapP(nonVerts);
sorted = sort([vmapMF,vmapPF],2);
[v i] = unique(sorted(:,1)); pairs = sorted(i,:);
for i = 1:size(pairs,1) 
    vmapT(pairs(i,:)) = Npv+i;
end
NfpT = length(v);

% interior dofs are independent, order at end
NipK = nnz(vmapT==0);
vmapT(vmapT==0) = Npv + NfpT + (1:NipK);

% restriction operator to assemble system
R = sparse(vmapT,1:Np*K,1);

% give back new boundary nodes
vmapBT = vmapT(vmapB);

