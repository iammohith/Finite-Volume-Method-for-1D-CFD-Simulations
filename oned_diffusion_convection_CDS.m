% Code for 1D steady-state diffusion and convection using FVM (CDS Scheme)

% Material properties and problem definition
rho = 1;           % Density (kg/m^3)
u = 0.1;          % Velocity (m/s)
K = 0.1;          % Diffusion coefficient (kg/m·s)
L = 1;            % Length of the 1D problem (m)
n = 5;            % Number of node points
dx = L / n;       % Distance between points (m)

% Boundary conditions
phi_a = 1;        % Transport property at point A
phi_b = 0;        % Transport property at point B

% Calculate flux and diffusion
F = rho * u;      % Flux
D = K / dx;       % Diffusion term

% Define nodal positions for FVM
fvm_x = dx/2:dx:(L - (dx/2)); % Nodal distance values
x = zeros(n + 2, 1);          % Matrix for x distances
x(2, 1) = dx / 2;             % First node point
x(end) = L;                   % Last element equals length L

% Populate nodal positions
for i = 3:n + 1
    x(i, 1) = x(i - 1, 1) + dx; % Array of nodal points
end

% Preallocate coefficient arrays
ap = zeros(n, 1);  % Coefficient for central nodes
aw = zeros(n, 1);  % Coefficient for west neighbors
ae = zeros(n, 1);  % Coefficient for east neighbors
su = zeros(n, 1);  % Source term coefficients
sp = zeros(n, 1);  % Pressure coefficient

% Calculate coefficients
for i = 2:n
    aw(i, 1) = D + (F / 2);            % Coefficient for west neighbors
end
for i = 1:n - 1
    ae(i, 1) = D - (F / 2);            % Coefficient for east neighbors
end
for i = 2:n - 1
    ap(i, 1) = (D + (F / 2)) + (D - (F / 2)); % Coefficient for central nodes
end

% Boundary conditions for coefficients
for i = 1:n
    if i == 1  % Node 1 (boundary condition at A)
        su(i, 1) = (2 * D + F) * phi_a;  % Source term at node 1
        sp(i, 1) = -(2 * D + F);          % Pressure coefficient at node 1
        ap(i, 1) = (D - (F / 2)) + (2 * D + F); % Coefficient for node 1
    elseif i == n  % Node n (boundary condition at B)
        su(i, 1) = (2 * D - F) * phi_b;  % Source term at node n
        sp(i, 1) = -(2 * D - F);          % Pressure coefficient at node n
        ap(i, 1) = (D + (F / 2)) + (2 * D - F); % Coefficient for node n
    end
end

% Solution using the Thomas Algorithm (TDMA)
S = zeros(n, n); % Coefficient matrix
X = zeros(n, 1); % Solution vector
B = zeros(n, 1); % Source vector
phi = zeros(n + 2, 1); % Solution for transport property
a = zeros(n - 1, 1); % Lower diagonal coefficients
b = zeros(n, 1);     % Main diagonal coefficients
c = zeros(n - 1, 1); % Upper diagonal coefficients
d = zeros(n, 1);     % Right-hand side vector

% Assign boundary values
phi(1) = phi_a; % Value at boundary A
phi(end) = phi_b; % Value at boundary B

% Construct the coefficient matrix S and source vector B
for i = 1:n
    if i == 1  % Node 1
        S(i, i) = (D - (F / 2)) + (2 * D + F); % Coefficient for node 1
        b(i, 1) = S(i, i);                    % Main diagonal at node 1
        B(i, 1) = su(i, 1);                   % Source term at node 1
        d(i, 1) = B(i, 1);                    % Right-hand side at node 1
    elseif i == n  % Node n
        S(i, i) = (D + (F / 2)) + (2 * D - F); % Coefficient for node n
        b(i, 1) = S(i, i);                    % Main diagonal at node n
        B(i, 1) = su(i, 1);                   % Source term at node n
        d(i, 1) = B(i, 1);                    % Right-hand side at node n
    end
end

% Fill in the coefficient matrix S
for i = 2:n - 1
    S(i, i) = (D + (F / 2)) + (D - (F / 2)); % Coefficient for nodes 2 to n-1
    b(i, 1) = S(i, i);                       % Main diagonal for nodes 2 to n-1
end
for i = 1:n - 1
    S(i, i + 1) = -(D - (F / 2));            % Coefficient for east neighbors
    c(i, 1) = S(i, i + 1);                   % Upper diagonal
    S(i + 1, i) = -(D + (F / 2));            % Coefficient for west neighbors
    a(i, 1) = S(i + 1, i);                   % Lower diagonal
end

% Thomas Algorithm implementation for TDMA
N = length(d); % Number of rows
c(1) = c(1) / b(1); % Division by zero risk
d(1) = d(1) / b(1); % Division by zero risk
for i = 2:N - 1
    temp = b(i) - a(i) * c(i - 1);
    c(i) = c(i) / temp;
    d(i) = (d(i) - a(i) * d(i - 1)) / temp;
end

% Back substitution to solve for X
d(N) = (d(N) - a(N - 1) * d(N - 1)) / (b(N) - a(N - 1) * c(N - 1));
X(N, 1) = d(N);
for i = N - 1:-1:1
    X(i, 1) = d(i, 1) - c(i, 1) * X(i + 1, 1);
end

% Assign solution to transport property array
for i = 1:n
    phi(i + 1, 1) = X(i, 1); % Fill phi array with the computed values
end

% Plotting results
figure(1)
plot(x, phi, '*')
title('Transport Property Plot for 1D Steady State Conduction & Convection')
xlabel('X [m]')
ylabel('Transport Property [\phi]')
