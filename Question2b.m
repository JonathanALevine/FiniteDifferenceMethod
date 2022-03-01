close all;  
clear; %intialization

set(0,'DefaultFigureWindowStyle','docked');

save_plots = 1;

% Potential and Length and Width of the rectangular region
W = 2;
L = 3/2*W;
V0 = 1;

dx = linspace(0.01, 1, 20);
dy = linspace(0.01, 1, 20);

for counter = 1:length(dx)
    counter
    nx = floor(L/dx(counter));
    ny = floor(W/dy(counter));

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

    V=G\B';

    for i=1:nx
        for j=1:ny
            N=j+(i-1)*ny;
            VMap(j,i)=V(N);
        end
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

    % Show plot of the electric field
    [Ex, Ey] = gradient(VMap);
    Jx = -cMap'.*Ex;
    Jy = -cMap'.*Ey;

    Ex = -Ex;
    Ey = -Ey;

    eFlowx = cMap .* Ex';
    eFlowy = cMap .* Ey';

    C0 = sum(eFlowx(1,:));
    Cnx = sum(eFlowx(nx, :));
    Curr(counter) = (C0 + Cnx) * 0.5;
end

figure('name', 'Current vs. Mesh Spacing');
plot(dx, Curr)
hold on 
plot(dx, Curr, 'bo');
hold off
xlabel('Mesh Spacing');
ylabel('Current (A)');
title('Current vs. Mesh Spacing');
ylim([0 1.5*max(Curr)])
xlim([min(dx) max(dx)])

if save_plots
    FN2 = 'Question 2b';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end
