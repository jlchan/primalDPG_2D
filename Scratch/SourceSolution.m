% Plot the solution to the Helmholtz equation with a given source
clear all;

Box_x = 3;
Scale = 0.5;
Box_y = Box_x/Scale;

Nx = 100;
Ny = Nx/Scale;

wavenumber = 10;

XX = linspace(-Box_x, Box_x, Nx);
YY = linspace(-Box_y, Box_y, Ny);
hx = XX(2) - XX(1);
hy = YY(2) - YY(1);
[X, Y] = meshgrid(XX, YY);

Source_size  = 0.5;
Source_shift = 2;
Source =  max(Source_size^2 - X.^2-(Y-Source_shift).^2, 0) + max(Source_size^2 - X.^2-(Y+Source_shift).^2, 0) ;

% plot the source
figure(1); clf; hold on; axis equal; axis off;
imagesc(Source);


% plot the solution to the Helmholtz equation
I = sqrt(-1);
Field = 0*X;

[m, n] = size(Source);
for i=1:m
    i
    for j=1:n
        
        if Source(i, j) ~= 0
            
            x0 = X(i, j);
            y0 = Y(i, j);
            
            % add the contribution from the current source
            Field = Field + (I/4)*besselh(0, 1, wavenumber*sqrt((X-x0).^2+(Y-y0).^2) + eps)*Source(i, j)*hx*hy;
        end
        
    end
end

figure(2); %clf; hold on; 
axis equal; axis off;
% imagesc(real(Field));
surf(real(Field))
shading interp
