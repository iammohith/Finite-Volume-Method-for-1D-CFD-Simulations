% Code for 1D Steady-State Heat Conduction without Heat Generation using FVM (Finite Volume Method)

% Material and problem properties
k = 1000; % Thermal conductivity of the material (W/m-K)
A = 10 * 10^-3; % Cross-sectional area perpendicular to heat flow (m^2)
l = 0.5; % Length of the 1D problem (m)
n = 5; % Number of nodes (grid points)
dx = l/n; % Distance between nodes (m)

% Boundary conditions
Ta = 100; % Temperature at the left boundary (A) (°C)
Tb = 500; % Temperature at the right boundary (B) (°C)

% Nodal distances for FVM
fvm_x = dx/2:dx:(l-(dx/2)); % Nodal points for FVM
x = zeros(n+2,1); % Preallocating array for nodal distances
x(2,1) = dx/2; % First node point
x(end) = l; % Last node is at the end of the length

% Constructing nodal distance array (x) for the grid
for i = 3:n+1
    x(i,1) = x(i-1,1) + dx; % Filling in the nodal distances
end

% Preallocating arrays for coefficients in FVM
ap = zeros(n,1); % Central coefficients
aw = zeros(n,1); % West (left) coefficients
ae = zeros(n,1); % East (right) coefficients
su = zeros(n,1); % Source term (SU)
sp = zeros(n,1); % Source coefficient (SP)

% Filling coefficients for nodes
for i = 2:n
    aw(i,1) = (k * A) / dx; % West coefficient (for nodes 2 to n)
end
for i = 1:n-1
    ae(i,1) = (k * A) / dx; % East coefficient (for nodes 1 to n-1)
end
for i = 2:n-1
    ap(i,1) = (2 * k * A) / dx; % Central coefficient (for nodes 2 to n-1)
end

% Handling boundary nodes (1 and n)
for i = 1:n
    if i == 1
        % At node 1 (left boundary)
        su(i,1) = ((2 * k * A) / dx) * Ta; % Source term
        sp(i,1) = -((2 * k * A) / dx); % Source coefficient
        ap(i,1) = (3 * k * A) / dx; % Central coefficient
    elseif i == n
        % At node n (right boundary)
        su(i,1) = ((2 * k * A) / dx) * Tb; % Source term
        sp(i,1) = -((2 * k * A) / dx); % Source coefficient
        ap(i,1) = (3 * k * A) / dx; % Central coefficient
    end
end

% TDMA (Tri-Diagonal Matrix Algorithm) solver setup
S = zeros(n,n); % Preallocating matrix for coefficients
X = zeros(n,1); % Preallocating solution vector for temperatures
B = zeros(n,1); % Preallocating source term vector
T = zeros(n+2,1); % Preallocating temperature vector with boundary values
a = zeros(n-1,1); % Lower diagonal for TDMA
b = zeros(n,1); % Main diagonal for TDMA
c = zeros(n-1,1); % Upper diagonal for TDMA
d = zeros(n,1); % Right-hand side for TDMA

% Setting boundary conditions in temperature array
T(1) = Ta; % Temperature at left boundary
T(end) = Tb; % Temperature at right boundary

% Constructing matrices for TDMA
for i = 1:n
    if i == 1
        % At node 1 (left boundary)
        S(i,i) = (3 * k * A) / dx; % Central coefficient
        b(i,1) = (3 * k * A) / dx;
        B(i,1) = ((2 * k * A) / dx) * Ta; % Source term
        d(i,1) = B(i,1);
    elseif i == n
        % At node n (right boundary)
        S(i,i) = (3 * k * A) / dx;
        b(i,1) = (3 * k * A) / dx;
        B(i,1) = ((2 * k * A) / dx) * Tb; % Source term
        d(i,1) = B(i,1);
    end
end

% For interior nodes (2 to n-1)
for i = 2:n-1
    S(i,i) = (2 * k * A) / dx; % Central coefficient for interior nodes
    b(i,1) = (2 * k * A) / dx;
end

% Filling in the coefficients for TDMA
for i = 1:n-1
    S(i,i+1) = -((k * A) / dx); % East coefficient
    c(i,1) = -((k * A) / dx);
    S(i+1,i) = -((k * A) / dx); % West coefficient
    a(i,1) = -((k * A) / dx);
end

% TDMA solver for system of equations
N = length(d); % Number of nodes

% Forward elimination phase
c(1) = c(1) / b(1); % Normalize first row
d(1) = d(1) / b(1); % Normalize RHS
for i = 2:N-1
    temp = b(i) - a(i) * c(i-1); % Modify the main diagonal
    c(i) = c(i) / temp; % Update upper diagonal
    d(i) = (d(i) - a(i) * d(i-1)) / temp; % Update RHS
end

% Backward substitution
d(N) = (d(N) - a(N-1) * d(N-1)) / (b(N) - a(N-1) * c(N-1)); % Solve for the last node
X(N,1) = d(N);
for i = N-1:-1:1
    X(i,1) = d(i,1) - c(i,1) * X(i+1,1); % Solve for remaining nodes
end

% Update temperature values with solved interior temperatures
for i = 1:n
    T(i+1,1) = X(i,1);
end

% Plotting the temperature distribution
figure (1)
plot(x,T,'*')
title('Temperature Distribution for 1D Steady-State Heat Conduction without Heat Generation')
xlabel('X [m]')
ylabel('Temperature [°C]')
