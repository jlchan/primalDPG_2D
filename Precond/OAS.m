function b = OAS(x,Af,Afi)

d = zeros(size(x));
for i = 1:length(Afi)
    d(Afi{i}) = d(Afi{i}) + 1;
end
d = 1./d; % scale 

b = zeros(size(x));
for i = 1:length(Afi)
    b(Afi{i}) = b(Afi{i}) +  d(Afi{i}).*(Af{i}\x(Afi{i}));
end

