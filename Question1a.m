close all;  
clear; %intialization

set(0,'DefaultFigureWindowStyle','docked')

save_plots = 1;

% Potential and Length and Width of the rectangular region
W = 2;
L = 3/2*W;
V0 = 1;

dx = 0.1;
dy = 0.1;

nx = L/dx;
ny = W/dy;

% GV = F form
%  The G matrix
G = sparse(nx*ny);
F = zeros(nx*ny, 1);

% The terms that fill in the G matrix are
term1 = -2*(1/dx^2+1/dy^2);
term2 = 1/dx^2;
term3 = 1/dy^2;

% Fill in the F vector
for y=1:ny
    n = MapNode(1, y, nx);
    F(n) = V0;
end

% Fill in the G matrix
for x=2:(nx-1)
    for y=2:(ny-1)
        n = MapNode(x,y,nx);
        G(n,n) = term1;
        G(n,MapNode(x-1,y,nx)) = term2;
        G(n,MapNode(x+1,y,nx)) = term2;
        G(n,MapNode(x,y-1,nx)) = term3;
        G(n,MapNode(x,y+1,nx)) = term3;
    end
end
for y=1:ny
    n = MapNode(1,y,nx);
    G(n,n) = 1;
    n = MapNode(nx,y,nx);
    G(n,n) = 1;
end
for x=2:(nx-1)
    n = MapNode(x,1,nx);
    G(n,n) = 1;
    G(n,MapNode(x,2,nx)) = -1;
    n = MapNode(x,ny,nx);
    G(n,n) = 1;
    G(n,MapNode(x,ny-1,nx)) = -1;
end

V = G\F;
V = reshape(V,[],ny)';

figure('name', 'Question 1(a) Solution');
surf(linspace(0,L,nx),linspace(0,W,ny),V);
xlabel('x');
ylabel('y');
zlabel('Potential (V)')
title(sprintf('FD Solution', dx));
set(gca, 'View', [45 45]);
c = colorbar;
c.Label.String = 'Potential Scale (V)';

if save_plots
    FN2 = 'Question 1a - FD Solution';   
    print(gcf, '-dpng', '-r600', FN2);  %Save graph in PNG
end

