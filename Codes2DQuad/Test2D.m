Globals2D;

for N=3:20

r = JacobiGL(0, 0, N)*ones(1,N+1);
s = r';
r = r(:); s = s(:);

x2 = .8; x3 = 1.4; y3 = 1.7;
y2 = -.7; x4 = -1.5; y4 = 1.2;
VX = [-1 x2 x3 x4];
VY = [-1 y2 y3 y4];

EToV = [1 2 3 4];
va = EToV(:,1)'; vb = EToV(:,2)'; 
vc = EToV(:,3)'; vd = EToV(:,4)'; 


x = 0.25*((1-r).*(1-s)*VX(va)+...
	   (1+r).*(1-s)*VX(vb)+...
	   (1+r).*(1+s)*VX(vc)+...
	   (1-r).*(1+s)*VX(vd));

y = 0.25*((1-r).*(1-s)*VY(va)+...
	   (1+r).*(1-s)*VY(vb)+...
	   (1+r).*(1+s)*VY(vc)+...
	   (1-r).*(1+s)*VY(vd));

figure(1)

plot(reshape(x, N+1,N+1), reshape(y, N+1, N+1))

hold on
plot(reshape(x, N+1,N+1)', reshape(y, N+1, N+1)')
scatter(x(:), y(:), 'r*')
hold off


figure(2)
% out of origin
r = JacobiGL(0, 0, N);

I = ones(1,N+1);
Z = zeros(1,N+1);
a1s = [ I; -I; Z];
b1s = [-I; r'; Z];

a2s = [ -I; I; Z];
b2s = [-r';-I; Z];

newr = [];
news = [];

for n=1:N+1
  for m=1:N+1
    
    [test, xint] = EdgeIntersect3D(a1s(:,n), b1s(:,n), a2s(:,m), b2s(:,m));
    
    if(test==true)
      newr = [newr;xint(1)];
      news = [news;xint(2)];
    end
  end
end

newr = [newr(2:end-1); r];
news = [news(2:end-1);-r];

newr = [2*rand(1000,1)-1];
news = [2*rand(1000,1)-1];

ids = find(newr+news<=0);
newr = newr(ids);
news = news(ids);

for loop=1:1000
  if(~mod(loop,1000))
    fprintf(1,'. ');
   end
    
  newr = [newr;2*rand(1000,1)-1];
  news = [news;2*rand(1000,1)-1];
  
  ids = find(newr+news<=0);
  newr = newr(ids);
  news = news(ids);

  [newrs] = unique([newr,news], 'rows');
  newr = newrs(:,1);
  news = newrs(:,2);
  
  % build Vandermonde2D
  V = Vandermonde2D(N, newr, news);
  
  [Dr, Ds] = Dmatrices2D(N, newr, news, V);
  
  ids = find(max(abs(Dr), [], 1)>1e-10);
  
  % choose best fits
  newr = newr(ids);
  news = news(ids);

end

V = Vandermonde2D(N, newr, news);

[Dr, Ds] = Dmatrices2D(N, newr, news, V);

Nsamples = 10000;
Nlevels = 10;
maxleb(N) = Lebesgue2D(N, newr, news, Nsamples, Nlevels)

VX = [-1, 1, 0];
VY = [ 0, 0, 2*cos(pi/3)];

r = newr;
s = news;
EToV = [1,2,3];
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

scatter(x(:), y(:), 'r*')
hold on
plot( [a1s(1,:);b1s(1,:)], [a1s(2,:);b1s(2,:)], 'k-');
plot( [a2s(1,:);b2s(1,:)], [a2s(2,:);b2s(2,:)], 'k-');
hold off

sk = 1;
for n=0:N
  for m=0:N-n
    fn = (newr.^n).*(news.^m);
    dfndr = n*(newr.^(n-1)).*(news.^m);
    dfnds = m*(newr.^(n)).*(news.^(m-1));
    err(sk,:) = [max(abs(dfndr-Dr*fn)),    max(abs(dfnds-Ds*fn))];
    sk = sk+1;

  end
end

end
