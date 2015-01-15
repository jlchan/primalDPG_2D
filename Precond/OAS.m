function b = OAS(x,Af)
nSubD = size(Af,1);

b = zeros(size(x));
for i = 1:nSubD
    %     noise = randn(size(Af{i,2},2));
    %     noise = 1e-6*noise'*noise;
    %     sol = (Af{i,2}+noise)\x(Af{i,1});
    
    sol = Af{i,2}\x(Af{i,1});
    
    %[L U] = luinc(Af{i,2},1e-4);
    %sol = U\(L\x(Af{i,1}));
%     sol = gmres(Af{i,2},x(Af{i,1}),size(Af{i,2},2),1e-4,25);
    
    b(Af{i,1}) = b(Af{i,1}) +  sol; % shouldn't scale overlapping contributions
end
