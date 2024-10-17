% Code for 1D steady-state heat conduction in a fin using the Finite Volume Method (FVM)

% Parameters
L = 1;                  % Length of the fin (m)
m = 25;                 % m = h*P/(K*A), where h = convective heat transfer coefficient
                         % P = perimeter, K = thermal conductivity, A = cross-sectional area
Tb = 100;               % Base temperature of the fin (°C)
T_inf = 20;             % Ambient temperature (°C)
n = 10;                 % Number of node points
dx = L/n;               % Distance between nodes (m)

% Initialize nodal distance values
fvm_x = dx/2:dx:(L-(dx/2)); % Nodal distances for FVM
x = zeros(n+1, 1);          % Array for nodal distances
x(2, 1) = dx/2;             % First node point used in FVM

% Fill in the x array with nodal distances
for i = 3:n+1
    x(i, 1) = x(i-1, 1) + dx; % Cumulative distance values
end

% Preallocate arrays for coefficients
ap = zeros(n, 1); % Coefficient for the main diagonal
aw = zeros(n, 1); % Coefficient for the west neighbor
ae = zeros(n, 1); % Coefficient for the east neighbor
su = zeros(n, 1); % Source term coefficient
sp = zeros(n, 1); % Coefficient for the right-hand side of the equation

% Calculate coefficients for the interior nodes
for i = 2:n
    aw(i, 1) = 1/dx;            % Coefficient for west neighbor
    ae(i, 1) = 1/dx;            % Coefficient for east neighbor
    ap(i, 1) = (2/dx) + (m*dx); % Coefficient for the main diagonal
    su(i, 1) = m*dx*T_inf;      % Source term
    sp(i, 1) = -m*dx;           % Coefficient for the right-hand side
end

% Adjust coefficients for the boundary nodes
% Node 1 (base of the fin)
ap(1, 1) = (3/dx) + (m*dx);               % Coefficient for the main diagonal
su(1, 1) = (m*dx*T_inf) + ((2/dx)*Tb);   % Source term
sp(1, 1) = -((m*dx) + (2/dx));            % Coefficient for the right-hand side

% Node n (tip of the fin)
ap(n, 1) = (1/dx) + (m*dx);               % Coefficient for the main diagonal

% Set up the system of equations using the TDMA (Tridiagonal Matrix Algorithm)
S = zeros(n, n);     % Matrix for coefficients
X = zeros(n, 1);     % Solution vector
B = zeros(n, 1);     % Right-hand side vector
T = zeros(n+1, 1);   % Temperature array
a = zeros(n-1, 1);   % Lower diagonal coefficients
b = zeros(n, 1);     % Main diagonal coefficients
c = zeros(n-1, 1);   % Upper diagonal coefficients
d = zeros(n, 1);     % Right-hand side vector for TDMA

T(1) = Tb;           % Set temperature at the base of the fin

% Fill in the coefficients for the first and last nodes
for i = 1:n
    if i == 1
        S(i, i) = ap(1, 1);          % Coefficient for node 1
        b(i, 1) = ap(1, 1);          % Main diagonal value for node 1
        B(i, 1) = su(1, 1);          % Source term for node 1
        d(i, 1) = su(1, 1);          % Right-hand side for node 1
    elseif i == n
        S(i, i) = ap(n, 1);          % Coefficient for node n
        b(i, 1) = ap(n, 1);          % Main diagonal value for node n
        B(i, 1) = su(n, 1);          % Source term for node n
        d(i, 1) = su(n, 1);          % Right-hand side for node n
    end
end

% Fill in the coefficients for the interior nodes
for i = 2:n-1
    S(i, i) = ap(i, 1);              % Coefficient for interior nodes
    b(i, 1) = ap(i, 1);              % Main diagonal value for interior nodes
    B(i, 1) = su(i, 1);              % Source term for interior nodes
    d(i, 1) = su(i, 1);              % Right-hand side for interior nodes
end

% Fill in the tridiagonal structure of the matrix
for i = 1:n-1
    S(i, i+1) = -ae(i, 1);           % Upper diagonal coefficients
    c(i, 1) = -ae(i, 1);             % Upper diagonal coefficients
    S(i+1, i) = -aw(i+1, 1);         % Lower diagonal coefficients
    a(i, 1) = -aw(i+1, 1);           % Lower diagonal coefficients
end

% Forward elimination for TDMA
N = length(d); % Number of rows
c(1) = c(1) / b(1); % Avoid division by zero
d(1) = d(1) / b(1); % Avoid division by zero
for i = 2:N-1
    temp = b(i) - a(i) * c(i-1);
    c(i) = c(i) / temp;
    d(i) = (d(i) - a(i) * d(i-1)) / temp;
end

% Back substitution for TDMA
d(N) = (d(N) - a(N-1) * d(N-1)) / (b(N) - a(N-1) * c(N-1));
X(N, 1) = d(N);
for i = N-1:-1:1
    X(i, 1) = d(i, 1) - c(i, 1) * X(i + 1, 1);
end

% Fill temperature array with solutions
for i = 1:n
    T(i + 1, 1) = X(i, 1); % Fill in temperatures at each node
end

% Plotting the results
figure(1)
plot(x, T, '*')
title('Temperature Distribution for 1D Steady-State Conduction in a Fin')
xlabel('Position along the fin (X) [m]')
ylabel('Temperature [°C]')
