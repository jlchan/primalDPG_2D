% plots solution and nodal values. requires u = size Np*K

function plotSol(u,Nplot)

Globals2D;

if nargin<2
    Nplot = 2*N;
end

va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
[xu,yu] = EquiNodes2D(Nplot); [ru, su] = xytors(xu,yu);
Vu = Vandermonde2D(N,ru,su); Iu = Vu*invV;
xu = 0.5*(-(ru+su)*VX(va)+(1+ru)*VX(vb)+(1+su)*VX(vc));
yu = 0.5*(-(ru+su)*VY(va)+(1+ru)*VY(vb)+(1+su)*VY(vc));
figure
color_line3(xu,yu,Iu*reshape(u,Np,K),Iu*reshape(u,Np,K),'.')
hold on
plot3(x(:),y(:),u(:),'o','markersize',8)