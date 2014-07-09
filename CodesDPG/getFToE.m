% gets face-to-element connections
% input: fpairs dimension: (unique faces x 2) - gives glob->unique faces
% output: FToE = unique faces x 2. FToE_i1,2 = elements that face i
% connects to

function FToE = getFToE(fpairs)

Globals2D;
NfacesU = size(fpairs,2);
FToE = zeros(NfacesU,2); % two elements on each side of a face
for k = 1:K
    for f = 1:Nfaces
        findex = (k-1)*Nfaces + f; 
        FToE(fpairs(1,:)==findex,1) = k;
        FToE(fpairs(2,:)==findex,2) = k;        
    end
end

