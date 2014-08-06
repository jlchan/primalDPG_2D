function plotFlux(uhat)
Globals2D;FaceGlobals2D;
uh = reshape(uhat,Nfrp,NfacesU);

figure
hold on
for i = 1:NfacesU
%     color_line3(xf(:,i),yf(:,i),uh(:,i),uh(:,i));
    plot3(xf(:,i),yf(:,i),uh(:,i));
    color_line3(xf(:,i),yf(:,i),uh(:,i),uh(:,i),'.');
end
