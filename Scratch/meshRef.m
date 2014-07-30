% mesh construction from [Nv, VX, VY, K, EToV]
Globals2D

[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell1.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('squarereg.neu');
N = 1;
StartUp2D;
% plotVerts

refLevel = zeros(K,1); % refinement level - start all cells off at 0
redGreen = ones(K,1); % start all cells as green

nRefs = 1;
ref = round(rand(nRefs,1)*K);
ref = [1; 8]
Nrefine = length(ref);
% nbrs = EToE(toRef,:);
% nbrs = setdiff(unique(nbrs(:)),toRef)
% nRef = length(toRef);
% for i = 1:nRef
%     nbrV = EToV(nbrs,:);
%     refV = EToV(toRef(i),:);    
% end
refineflag = zeros(K,1);
refineflag(ref) = 1;

%% nonconforming refinements!

newVIds = HCrefine2D(refineflag);

%% mark greens/close mesh

nbrs = EToE(ref,:);

green = zeros(K,1);
for i = 1:Nrefine            
    for j = 1:3
        % bisection
        if nbrs(i,j)~=ref(i)  % if neighbor is not refined node...  
            k = nbrs(i,j);
            newVert = newVIds(i,j);
            ncFace = EToF(ref(i),j);
            ncVIds = [ncFace mod(ncFace,Nfaces)+1];
            oppVId = setdiff(1:3,ncVIds);
            ncV = EToV(k,ncVIds);
            oppV = EToV(k,oppVId);
            T1 = [oppV ncV(1) newVert];
            T2 = [newVert ncV(2) oppV];
            if (green(k)==0) % if we can 
                EToV(k,:) = T1; % replace original triangle
                EToV(K+1,:) = T2;
                green(k) = k; % mark
                green(K+1) = k; % mark as bisected
                K = K+1;
            else
                toMerge = find(green==green(k));
                vMerge = [EToV(toMerge(1),1:2) EToV(toMerge(2),2)];
                EToV(toMerge(1),:) = vMerge;
                EToV(toMerge(2),:) = [];                                
                K = K-1;
                keyboard
                green(k)=0;
%                 refineflag = zeros(K,1); refineflag(toMerge(1)) = 1;
%                 HCrefine2D(refineflag);
                % remove EToV rows, extra vertices, and use hrefine to
                % close.
                
            end
            clf;plotVerts
%             keyboard
        end
    end
end


%% check


