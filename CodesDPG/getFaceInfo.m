% fM/fP = get unique node mapping - volume nodes to face nodes
% fpairs = 2 x NfacesU matrix of face numberings and corresponding matching
% face (in the second row of the matrix)

function [fM fP fpairs] = getFaceInfo()

Globals2D

fM = reshape(vmapM,Nfp,Nfaces*K); fP = reshape(vmapP,Nfp,Nfaces*K); 

% sort by node id in face and then get unique pairs of columns 
[tf, loc] = ismember(sort(fM,1)',sort(fP,1)','rows'); 
fpairs = [(1:length(loc))' loc(:)];
fpairs = unique(sort(fpairs,2),'rows')'; 

% fpairs(2,:) are duplicated face nodes, ignore them
fM = fM(:,fpairs(1,:)); fP = fP(:,fpairs(1,:));
