% Code for 1D Steady-State Diffusion and Convection using FVM (CDS Scheme)

% Physical parameters
rho = 1;               % Density (kg/m^3)
u = 0.1;               % Velocity (m/s)
K = 0.1;               % Diffusion coefficient (kg/m.s)
L = 1;                 % Length of the 1D domain (m)
n = 5;                 % Number of node points

% Discretization parameters
dx = L / n;            % Distance between points (m)
phi_a = 1;             % Transport property at point A
phi_b = 0;             % Transport property at point B
F = rho * u;           % Flux (kg/(m^2.s))
D = K / dx;            % Diffusion (m^2/s)

% Nodal distance values used for FVM
fvm_x = dx/2:dx:(L - (dx/2));
x = zeros(n + 2, 1);   % Initialize x distance matrix
x(2) = dx / 2;         % First node point used in FVM
x(end) = L;            % Last element is the length of the domain

% Generate nodal distances array
for i = 3:n + 1
    x(i) = x(i - 1) + dx; % Create array of nodal point distances
end

% Preallocate arrays for coefficients
ap = zeros(n, 1);      % Coefficients for pressure
aw = zeros(n, 1);      % Coefficients for west (left)
ae = zeros(n, 1);      % Coefficients for east (right)
su = zeros(n, 1);      % Source terms
sp = zeros(n, 1);      % Source terms for pressure

% Determine coefficients based on flow direction
if u > 0 
    % For positive flow
    for i = 2:n
        aw(i) = D + F;  % Coefficient for west
    end
    for i = 1:n - 1
        ae(i) = D;      % Coefficient for east
    end
    for i = 2:n - 1
        ap(i) = (D + F) + D;  % Coefficient for pressure
    end
    for i = 1:n
        if i == 1
            su(i) = (2 * D + F) * phi_a;  % Source term at node 1
            sp(i) = -(2 * D + F);          % Source term for pressure at node 1
            ap(i) = (3 * D + F);           % Coefficient for pressure at node 1
        elseif i == n
            su(i) = (2 * D) * phi_b;      % Source term at node n
            sp(i) = -(2 * D);              % Source term for pressure at node n
            ap(i) = (3 * D + F);           % Coefficient for pressure at node n
        end
    end
elseif u < 0
    % For negative flow
    for i = 2:n
        aw(i) = D;  % Coefficient for west
    end
    for i = 1:n - 1
        ae(i) = D - F;  % Coefficient for east
    end
    for i = 2:n - 1
        ap(i) = (2 * D - F);  % Coefficient for pressure
    end
    for i = 1:n
        if i == 1
            su(i) = (2 * D) * phi_a;       % Source term at node 1
            sp(i) = -(2 * D);               % Source term for pressure at node 1
            ap(i) = (3 * D - F);            % Coefficient for pressure at node 1
        elseif i == n
            su(i) = (2 * D - F) * phi_b;   % Source term at node n
            sp(i) = -(2 * D - F);           % Source term for pressure at node n
            ap(i) = (3 * D - F);            % Coefficient for pressure at node n
        end
    end
end

% Solution using the Thomas Algorithm (TDMA)
S = zeros(n, n);       % Coefficient matrix for the system
X = zeros(n, 1);       % Solution vector
B = zeros(n, 1);       % Right-hand side vector
phi = zeros(n + 2, 1); % Array to store transport property values

% Initialize boundary conditions
phi(1) = phi_a;        % Boundary condition at A
phi(end) = phi_b;      % Boundary condition at B

% Construct the system of equations based on flow direction
if u > 0
    for i = 1:n
        if i == 1
            S(i, i) = (3 * D + F);  % Coefficient for node 1
            B(i) = (2 * D + F) * phi_a;  % Source term at node 1
        elseif i == n
            S(i, i) = (3 * D + F);  % Coefficient for node n
            B(i) = (2 * D) * phi_b;  % Source term at node n
        end
    end
    for i = 2:n - 1
        S(i, i) = (2 * D + F);  % Coefficient for nodes 2 to n-1
    end
    for i = 1:n - 1
        S(i, i + 1) = -D;       % Coefficient for east
        S(i + 1, i) = -(D + F); % Coefficient for west
    end
elseif u < 0
    for i = 1:n
        if i == 1
            S(i, i) = (3 * D - F);  % Coefficient for node 1
            B(i) = (2 * D) * phi_a;  % Source term at node 1
        elseif i == n
            S(i, i) = (3 * D - F);  % Coefficient for node n
            B(i) = (2 * D - F) * phi_b;  % Source term at node n
        end
    end
    for i = 2:n - 1
        S(i, i) = (2 * D - F);  % Coefficient for nodes 2 to n-1
    end
    for i = 1:n - 1
        S(i, i + 1) = -(D - F);  % Coefficient for east
        S(i + 1, i) = -D;        % Coefficient for west
    end
end

% Prepare for TDMA
N = length(B); % Number of equations
c(1) = c(1) / b(1); % First entry adjustment for TDMA
d(1) = d(1) / b(1); % First entry adjustment for TDMA

% Forward elimination
for i = 2:N - 1
    temp = b(i) - a(i) * c(i - 1);
    c(i) = c(i) / temp;
    d(i) = (d(i) - a(i) * d(i - 1)) / temp;
end

% Back substitution
d(N) = (d(N) - a(N - 1) * d(N - 1)) / (b(N) - a(N - 1) * c(N - 1));
X(N) = d(N);
for i = N - 1:-1:1
    X(i) = d(i) - c(i) * X(i + 1);
end

% Store solution
phi(2:n + 1) = X;

% Plotting the results
figure(1);
plot(x, phi, '*');
title('Transport Property Plot for 1D Steady-State Conduction & Convection');
xlabel('X [m]');
ylabel('Transport Property [\phi]');
