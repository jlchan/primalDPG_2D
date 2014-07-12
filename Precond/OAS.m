function b = OAS(x,Af)
nSubD = size(Af,1);
% d = zeros(size(x));
% for i = 1:nSubD
%     d(Af{i,1}) = d(Af{i,1}) + 1;
% end
% d = 1./d; % scale
% d = d.^0; % don't scale

b = zeros(size(x));
for i = 1:nSubD
%         b(Af{i,1}) = b(Af{i,1}) +  d(Af{i,1}).*(Af{i,2}\x(Af{i,1}));
    b(Af{i,1}) = b(Af{i,1}) +  Af{i,2}\x(Af{i,1}); % don't scale - shouldn't scale overlapping contributions!
end
