close all;  
clear; %intialization

set(0,'DefaultFigureWindowStyle','docked')

nx = 100;
ny = 100;
iterations = 10000;

V = zeros(nx, ny);
V(:,1) = 1;
V(:, ny) = 0;

ForceBoundaryToZero = 0;

% Show plot of potential
figure(1)
for iteration = 1:iterations
    if ForceBoundaryToZero
        for x = 2:nx-1
            for y = 2:ny-1
                V(x, y) = (V(x+1, y) + V(x-1, y) + V(x, y+1) + V(x, y-1))/4;
            end
        end
    else
        for x = 2:nx-1
            V(:,x) = (V(:,x+1) + V(:,x-1))/2;
        end
    end

    if mod(iteration, 50) == 0
        image = imboxfilt(V, 3);
        imagesc(image);
        pause(0.05);
    end
end

% Show plot of the electric field
[Ex, Ey] = gradient(V);

figure(2)
quiver(-Ex', -Ey', 1); 