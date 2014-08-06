% refines elements using recursive bisection. 
% stolen from Tim Warburton's codes. Probably still some cruft, but works
% for triangle meshes.  

function Refine2D(eflag)

Globals2D;FaceGlobals2D;

% find face centers
centr = [0;0;-1];
cents = [-1;0;0];
Icent = Vandermonde2D(N, centr, cents)/Vandermonde2D(N, r, s);
centx = (Icent*x)';
centy = (Icent*y)';

% base rotation
baserot = [ 1 2 3; 2 3 1; 3 1 2];

% mark bases
for k1=1:K
    % find base
    b = FindBase2D(k1);
    
    % find rotation
    br = baserot(b,:);
    
%     % 2. rotate boundary types
%     BCType(k1,:) = BCType(k1,br);
    
    % 3. rotate vertices
    EToV(k1,:) = EToV(k1,br);
    
    % 4. rotate flags
    eflag(k1,:) = eflag(k1,br);
    
    % 5. rotate face centers
    centx(k1,:) = centx(k1,br);
    centy(k1,:) = centy(k1,br);
    
end

[EToE,EToF] = tiConnect2D(EToV);

% make sure refinement is continuous
for k1=1:K
    for f1=1:Nfaces
        k2 = EToE(k1,f1); f2 = EToF(k1,f1);
        eflag(k1,f1) = max(eflag(k1,f1), eflag(k2,f2));
    end
end

% now keep propagating refinements of elements
% until all elements to be refined have marked bases
while(1)
    cnt = 0;
    for k1=1:K
        if(eflag(k1,1)==0)
            if(sum(eflag(k1,:))~=0)
                eflag(k1,1) = 1;
                k2 = EToE(k1,1);
                f2 = EToF(k1,1);
                eflag(k1,1)  = 1;
                eflag(k2,f2) = 1;
                cnt = cnt+1;
            end
        end
    end
    
    if(cnt==0)
        break;
    end
end

% label mid edges
EToM = reshape(1:Nfaces*K, Nfaces, K)';

for k1=1:K
    for f1=1:Nfaces
        k2 = EToE(k1,f1);
        f2 = EToF(k1,f1);
        EToM(k1,f1) = max(EToM(k1,f1),EToM(k2,f2));
        EToM(k2,f2) = EToM(k1,f1);
    end
end
EToM = EToM + length(VX) + 1;

vnum = [ 1 2; 2 3 ; 3 1 ];
newVX = VX;
newVY = VY;
newVX(EToM) = centx;
newVY(EToM) = centy;
if(0)
    for k1=1:K
        for f1=1:Nfaces
            if(0)
                x1 = VX(EToV(k1,vnum(f1,1)));
                y1 = VY(EToV(k1,vnum(f1,1)));
                x2 = VX(EToV(k1,vnum(f1,2)));
                y2 = VY(EToV(k1,vnum(f1,2)));
                
                newVX(EToM(k1,f1)) = 0.5*(x1+x2);
                newVY(EToM(k1,f1)) = 0.5*(y1+y2);
            end
            
            newVX(EToM(k1,f1)) = centx(k1,f1);
            newVY(EToM(k1,f1)) = centy(k1,f1);
        end
    end
end

% now refine
newK = sum(sum(eflag))+K;
newEToE   = zeros(newK, Nfaces);
newEToV   = zeros(newK, Nfaces);
% newBCType = zeros(newK, Nfaces);

sk = 1;
for k1=1:K
    sflag = sum(eflag(k1,:));
    
    v1 =   EToV(k1,1); v2 =   EToV(k1,2); v3 =   EToV(k1,3);
    m1 =   EToM(k1,1); m2 =   EToM(k1,2); m3 =   EToM(k1,3);
%     b1 = BCType(k1,1); b2 = BCType(k1,2); b3 = BCType(k1,3);
    
    if(sflag==0)
        newEToV(sk,:)   = [v1,v2,v3];
%         newBCType(sk,:) = [b1,b2,b3];
        sk = sk+1;
    elseif(sflag==1)
        newEToV(sk,:)   = [v1,m1,v3];
%         newBCType(sk,:) = [b1,0,b3];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v2,v3];
%         newBCType(sk,:) = [b1,b2,0];
        sk = sk+1;
    elseif(sflag==3)
        newEToV(sk,:)   = [v1,m1,m3];
%         newBCType(sk,:) = [b1,0,b3];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v2,m2];
%         newBCType(sk,:) = [b1,b2,0];
        sk = sk+1;
        newEToV(sk,:)   = [m1,m2,v3];
%         newBCType(sk,:) = [0,b2,0];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v3,m3];
%         newBCType(sk,:) = [0,b3,0];
        sk = sk+1;
    elseif(eflag(k1,2)==1)
        newEToV(sk,:)   = [v1,m1,v3];
%         newBCType(sk,:) = [b1,0,b3];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v2,m2];
%         newBCType(sk,:) = [b1, b2, 0];
        sk = sk+1;
        newEToV(sk,:)   = [m1,m2,v3];
%         newBCType(sk,:) = [0,b2,0];
        sk = sk+1;
    else
        newEToV(sk,:)   = [v1,m1,m3];
%         newBCType(sk,:) = [b1,0,b3];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v3,m3];
%         newBCType(sk,:) = [0,b3,0];
        sk = sk+1;
        newEToV(sk,:)   = [m1,v2,v3];
%         newBCType(sk,:) = [b1,b2,0];
        sk = sk+1;
    end
end

K      = newK;
EToV   = newEToV;
% BCType = newBCType;
[EToE,EToF] = tiConnect2D(EToV);
VX = newVX;
VY = newVY;

% relabel VX, VY
[uV i j] = unique(EToV);
nV = length(uV);
% map = [(1:nuV)' uV];

for i = 1:nV
    EToV(EToV==uV(i)) = i;
end
VX = VX(uV); VY = VY(uV);

StartUp2D; FaceStartUp2D;

%  PlotMesh2D();
%   ha= trimesh(EToV, VX, VY);
%   set(ha, 'Color', 'black')
%   axis equal

function b = FindBase2D(k)

Globals2D;

x1 = VX(EToV(k,1)); y1 = VY(EToV(k,1));
x2 = VX(EToV(k,2)); y2 = VY(EToV(k,2));
x3 = VX(EToV(k,3)); y3 = VY(EToV(k,3));

L2 = [(x1-x2).^2+(y1-y2).^2, (x2-x3).^2+(y2-y3).^2, (x3-x1).^2+(y3-y1).^2];

[maxL2,b] = max(L2);

