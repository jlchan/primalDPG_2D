function plotFlux(uhat)
Globals2D;FaceGlobals2D;

uh = reshape(uhat,Nfrp,NfacesU);
hold on
for i = 1:NfacesU
    color_line3(xf(:,i),yf(:,i),uh(:,i),uh(:,i));
    color_line3(xf(:,i),yf(:,i),uh(:,i),uh(:,i),'.');
end
