close all;  
clear; %intialization

set(0,'DefaultFigureWindowStyle','docked');

save_plots = 0;

% Potential and Length and Width of the rectangular region
W = 2;
L = 3/2*W;
V0 = 1;

dx = 0.05;
dy = 0.05;

nx = L/dx;
ny = W/dy;

% Conductivities
sigma_inside = 0.01;
sigma_outside = 1;

% GV = 0 form
G = sparse(nx*ny);
cMap = zeros(nx, ny);
B = zeros(1, nx*ny);

% Define the limits of the boxes
boxes = [nx*2/5 nx*3/5 ny*1/3 ny*2/3]; 

for i=1:nx
    for j=1:ny
        if i > boxes(1) && i < boxes(2) & (j < boxes(3) || j > boxes(4))
            cMap(i, j) = sigma_inside;
        else
            cMap(i, j) = sigma_outside;
        end
    end
end

for x = 1:nx
    for y = 1:ny
        n = MapNode(y, x, ny);
        nxp = MapNode(y, x+1, ny);
        nxm = MapNode(y, x-1, ny);
        nyp = MapNode(y+1, x, ny);
        nym = MapNode(y-1, x, ny);
        
        if x == 1
            G(n, n) = 1;
            B(n) = V0;

        elseif x == nx
            G(n, n) = 1;
            B(n) = 0;
            
        elseif y == 1
            G(n, nxp) = (cMap(x+1, y) + cMap(x,y))/2;
            G(n, nxm) = (cMap(x-1, y) + cMap(x,y))/2;
            G(n, nyp) = (cMap(x, y+1) + cMap(x,y))/2;            
            G(n, n) = -(G(n,nxp)+G(n,nxm)+G(n,nyp));
            
        elseif y == ny
            G(n, nxp) = (cMap(x+1, y) + cMap(x,y))/2;
            G(n, nxm) = (cMap(x-1, y) + cMap(x,y))/2;
            G(n, nym) = (cMap(x, y-1) + cMap(x,y))/2;
            G(n, n) = -(G(n,nxp)+G(n,nxm)+G(n,nym));
            
        else
            G(n, nxp) = (cMap(x+1, y) + cMap(x,y))/2;
            G(n, nxm) = (cMap(x-1, y) + cMap(x,y))/2;
            G(n, nyp) = (cMap(x, y+1) + cMap(x,y))/2;
            G(n, nym) = (cMap(x, y-1) + cMap(x,y))/2;
            G(n, n) = -(G(n,nxp)+G(n,nxm)+G(n,nyp)+G(n,nym));
        end
    end
end

% Sigma(x,y) Surface Plot
figure('name', 'sigma(x, y)')
surf(cMap');
xlabel("X position")
ylabel("Y position")
zlabel("\sigma(x, y)")
title("Conductivity")
axis tight
view(0,90)
c = colorbar;
c.Label.String = 'Conductivity Scale (S/m)';

if save_plots
    FN2 = 'Question 2a - sigma';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end

% Get the solution of the voltage
V = G\B';

VMap = zeros(ny, nx);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        VMap(j,i) = V(n);
    end
end

% V(x,y) Surface Plot
figure('name', 'V(x, y)')
surf(VMap)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")
title("Potential")
view(0, 90)
c = colorbar;
c.Label.String = 'Potential Scale (V)';

if save_plots
    FN2 = 'Question 2a - Potential';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end

% Show plot of the electric field
[Ex, Ey] = gradient(VMap);
Jx = -cMap'.*Ex;
Jy = -cMap'.*Ey;

figure('name', 'E(x, y)')
quiver(-Ex, -Ey)
xlabel("X position")
ylabel("Y position")
zlabel("E(x, y)")
title("Electric Field")
axis tight
view(0, 90);

if save_plots
    FN2 = 'Question 2a - EField';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end

figure('name', 'J(x, y)')
quiver(Jx, Jy)
xlabel("X position")
ylabel("Y position")
zlabel("J(x, y)")
title("Current Density")
axis tight
view(0, 90);


if save_plots
    FN2 = 'Question 2a - CurrentDensity';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap .* Ex';
eFlowy = cMap .* Ey';

C0 = sum(eFlowx(1,:));
Cnx = sum(eFlowx(nx, :));
Curr = (C0 + Cnx) * 0.5;

