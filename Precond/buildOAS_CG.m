function Mhandle = buildOAS_CG(R,A,Norder)

Globals2D;

if nargin<3
    Norder = N;
end

% build OAS preconditioner algebraically:
if (size(EToV,2)==4) % if quad
    Norderp = (Norder+1)^2;
else
   Norderp = (Norder+1)*(Norder+2)/2;
end

[rr cc] = find(R); 
rr = reshape(rr,Norderp,K);cc = reshape(cc,Norderp,K);

patches = 0;
if patches
    % build patches around vertices
    vertices = rr([1 Norder+1 Norderp],:);
    vnodes = unique(vertices(:));
    Nverts = length(vnodes);    
    Ak = cell(Nverts,2);
    for i = 1:length(vnodes)
        [~, elems] = find(vnodes(i)==vertices);
        inds = unique(rr(:,elems));        
        Ak{i,1} = inds;
        Ak{i,2} = A(inds,inds);          
    end
else
    % no overlap/face overlap only
    Ak = cell(K,2);
    for k = 1:K
        %   elems = unique([k EToE(k,:)]); % face overlap only
        elems = k; % no overlap
        inds = unique(rr(:,elems));
        Ak{k,1} = inds;
        Ak{k,2} = A(inds,inds);
%         if (k==30)
%             clf;
%             ti = cc(:,elems);
%             plotVerts;hold on
%             plot(x(ti),y(ti),'ro')
%             keyboard
%         end
    end
end

% coarse grid solver
[Rc1 Ir1 vmapBT1 xr1 yr1] = pRestrictCG(Norder,1); % interp down
Rp1 = Rc1*Ir1';
Rs = spdiag(1./sum(R,2))*R; % divide by shared nodal contributions - undo assembly
R1 = Rs*Rp1'; % interp down to P1
P1 = @(b) R1*((R1'*A*R1)\(R1'*b)); % should be able to "ignore" BC imposition - pcg acts on residual e(vmapBT) = 0

Mhandle = @(x) P1(x) + OAS(x,Ak); % coarse solver + OAS solver