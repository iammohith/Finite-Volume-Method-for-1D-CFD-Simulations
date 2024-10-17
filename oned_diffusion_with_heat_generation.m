% Code for 1D Steady State Heat Conduction with Heat Generation using FVM
% Material and geometric properties
k = 0.5;                   % Thermal conductivity of material (W/m-K)
A = 1;                     % Perpendicular area (m^2)
L = 0.02;                  % Length of 1D problem (m)
q = 1000 * 10^3;          % Heat generation rate (W/m^3)
n = 5;                     % Number of node points
dx = L / n;                % Distance between node points (m)

% Boundary temperatures
Ta = 100;                  % Temperature at A (°C)
Tb = 200;                  % Temperature at B (°C)

% Nodal distance values used for FVM
fvm_x = dx/2:dx:(L - (dx/2));
x = linspace(0, L, n + 2)';  % Creating a vector for nodal distances

% Preallocating coefficient arrays
ap = zeros(n, 1);           % Coefficient for internal nodes
aw = zeros(n, 1);           % Coefficient for west neighbor
ae = zeros(n, 1);           % Coefficient for east neighbor
su = zeros(n, 1);           % Source term
sp = zeros(n, 1);           % Source term for boundary conditions

% Coefficients for internal nodes
for i = 2:n
    aw(i) = (k * A) / dx;   % West coefficient
    ae(i - 1) = (k * A) / dx; % East coefficient
end

% Coefficients for internal nodes (ap and su)
for i = 2:n-1
    ap(i) = (2 * k * A) / dx; % Coefficient for internal nodes
    su(i) = q * A * dx;       % Source term for internal nodes
end

% Boundary conditions
for i = 1:n
    if i == 1
        su(i) = q * A * dx + ((2 * k * A) / dx) * Ta; % Source term at node 1
        sp(i) = -((2 * k * A) / dx);                   % Source term for node 1
        ap(i) = (3 * k * A) / dx;                      % Coefficient for node 1
    elseif i == n
        su(i) = q * A * dx + ((2 * k * A) / dx) * Tb; % Source term at node n
        sp(i) = -((2 * k * A) / dx);                   % Source term for node n
        ap(i) = (3 * k * A) / dx;                      % Coefficient for node n
    end
end

% Solution using TDMA (Tri-Diagonal Matrix Algorithm)
S = zeros(n, n);           % Preallocating matrix for system of equations
X = zeros(n, 1);           % Preallocating solution vector
B = zeros(n, 1);           % Preallocating right-hand side vector
T = zeros(n + 2, 1);       % Preallocating temperature vector
a = zeros(n - 1, 1);       % Sub-diagonal coefficients for TDMA
b = zeros(n, 1);           % Main diagonal coefficients for TDMA
c = zeros(n - 1, 1);       % Super-diagonal coefficients for TDMA
d = zeros(n, 1);           % Right-hand side vector for TDMA

% Assign boundary temperatures
T(1) = Ta;                  % Temperature at boundary A
T(end) = Tb;                % Temperature at boundary B

% Constructing the system of equations
for i = 1:n
    if i == 1
        S(i, i) = (3 * k * A) / dx; % Main diagonal at node 1
        b(i) = (3 * k * A) / dx;    % Main diagonal value at node 1
        B(i) = q * A * dx + ((2 * k * A) / dx) * Ta; % Right-hand side at node 1
        d(i) = B(i);                % Assign to d vector
    elseif i == n
        S(i, i) = (3 * k * A) / dx; % Main diagonal at node n
        b(i) = (3 * k * A) / dx;    % Main diagonal value at node n
        B(i) = q * A * dx + ((2 * k * A) / dx) * Tb; % Right-hand side at node n
        d(i) = B(i);                % Assign to d vector
    end
end

for i = 2:n-1
    S(i, i) = (2 * k * A) / dx;   % Main diagonal for internal nodes
    b(i) = (2 * k * A) / dx;      % Main diagonal value for internal nodes
    B(i) = q * A * dx;             % Right-hand side for internal nodes
    d(i) = B(i);                   % Assign to d vector
end

for i = 1:n-1
    S(i, i + 1) = -((k * A) / dx); % Coefficient for east neighbor
    c(i) = -((k * A) / dx);        % Super-diagonal coefficient
    S(i + 1, i) = -((k * A) / dx); % Coefficient for west neighbor
    a(i) = -((k * A) / dx);        % Sub-diagonal coefficient
end

% Forward elimination for TDMA
N = length(d);                % Total number of equations
c(1) = c(1) / b(1);           % First super-diagonal coefficient
d(1) = d(1) / b(1);           % First right-hand side value

for i = 2:N-1
    temp = b(i) - a(i) * c(i-1); % Temporary variable for elimination
    c(i) = c(i) / temp;          % Update super-diagonal coefficient
    d(i) = (d(i) - a(i) * d(i-1)) / temp; % Update right-hand side value
end

% Back substitution for TDMA
d(N) = (d(N) - a(N-1) * d(N-1)) / (b(N) - a(N-1) * c(N-1));
X(N) = d(N);                  % Last solution value

for i = N-1:-1:1
    X(i) = d(i) - c(i) * X(i + 1); % Back substitution
end

% Storing the temperature results in T vector
for i = 1:n
    T(i + 1) = X(i);              % Assign temperature values to T vector
end

% Plotting the results
figure(1)
plot(x, T, '*')
title('Temperature Plot for 1D Steady State Heat Conduction with Heat Generation')
xlabel('X [m]')
ylabel('Temperature [°C]')
