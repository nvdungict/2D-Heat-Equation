% Explicit Scheme 2D Heat Equation

clear;
clc;

% General input
TOL = 1e-6;      % Tolerance for convergence
temp = 100;      % Temperature boundary condition

% THIS SECTION CALCULATES TEMPs FOR N=11
n = 11;  % grid size
x = linspace(0, 1, n);
dx = x(2) - x(1);
y = x;
dy = dx;

T11 = zeros(n);
T11(1, 1:n) = 0;    % BOTTOM
T11(n, 1:n) = temp; % TOP
T11(1:n, 1) = 0;    % LEFT
T11(1:n, n) = 0;    % RIGHT

dt = dx^2 / 4;  % Time step
error = 1;  % Initial error
k = 0;  % Iteration counter

while error > TOL
    k = k + 1;
    Told = T11;
    for i = 2:n-1
        for j = 2:n-1
            T11(i, j) = dt * ((Told(i+1, j) - 2*Told(i, j) + Told(i-1, j)) / dx^2 ...
                            + (Told(i, j+1) - 2*Told(i, j) + Told(i, j-1)) / dy^2) ...
                            + Told(i, j);
        end
    end
    error = max(max(abs(Told - T11)));
end

figure;
pcolor(x, y, T11), shading interp, xlabel('x'), ylabel('y'), colorbar;
title('Temperature Distribution for 11x11 Grid (Explicit Scheme)');

% THIS SECTION CALCULATES TEMPs FOR N=21
n = 21;  % grid size
x = linspace(0, 1, n);
dx = x(2) - x(1);
y = x;
dy = dx;

T21 = zeros(n);
T21(1, 1:n) = 0;    % BOTTOM
T21(n, 1:n) = temp; % TOP
T21(1:n, 1) = 0;    % LEFT
T21(1:n, n) = 0;    % RIGHT

dt = dx^2 / 4;  % Time step
error = 1;  % Initial error
k = 0;  % Iteration counter

while error > TOL
    k = k + 1;
    Told = T21;
    for i = 2:n-1
        for j = 2:n-1
            T21(i, j) = dt * ((Told(i+1, j) - 2*Told(i, j) + Told(i-1, j)) / dx^2 ...
                            + (Told(i, j+1) - 2*Told(i, j) + Told(i, j-1)) / dy^2) ...
                            + Told(i, j);
        end
    end
    error = max(max(abs(Told - T21)));
end

figure;
pcolor(x, y, T21), shading interp, xlabel('x'), ylabel('y'), colorbar;
title('Temperature Distribution for 21x21 Grid (Explicit Scheme)');
