function b = OAS(x,Af)
nSubD = size(Af,1);

b = zeros(size(x));
for i = 1:nSubD
    b(Af{i,1}) = b(Af{i,1}) +  Af{i,2}\x(Af{i,1}); % shouldn't scale overlapping contributions
end
