% Code for 1D SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) Algorithm
% This implementation uses a staggered grid approach for solving fluid flow

% Grid and fluid properties
n = 4;          % Number of grid points
A = 1;          % Cross-sectional area (m^2)
rho = 1;        % Density (kg/m^3)
d_co = 1;       % Ratio of coefficients A(i,J)/a(i,J)

% Velocity and pressure arrays
u = zeros(n,1);         % Preallocating matrix for velocities (m/s)
p_cor = zeros(n,1);     % Preallocating matrix for pressure correction (Pa)
u(1,1) = 10;            % Inlet velocity (m/s)
p_cor(end) = 0;         % Outlet pressure condition (Gauge pressure in Pa)

% Initial guess for velocities
u_guess = [8; 11; 7];   % Guess values of u (m/s)

% Preallocating arrays for solving pressure correction equation using TDMA
S = zeros(n-1,n-1);     % Coefficient matrix
X = zeros(n-1,1);       % Solution vector for pressure correction
B = zeros(n-1,1);       % RHS vector for TDMA
a = zeros(n-2,1);       % Lower diagonal (TDMA)
b = zeros(n-1,1);       % Main diagonal (TDMA)
c = zeros(n-2,1);       % Upper diagonal (TDMA)
d = zeros(n-1,1);       % Modified RHS vector for TDMA

% Coefficients for discretized equations
aw = rho * d_co * A;    % West coefficient (aw)
ae = rho * d_co * A;    % East coefficient (ae)
ap = aw + ae;           % Central coefficient (ap)

% Constructing the matrix for pressure correction
for I = 1:n-1
    if I == 1
        % Boundary condition at inlet
        S(I,I) = ae;
        b(I,1) = ae;
        B(I,1) = -(rho * u_guess(I,1) * A) + (rho * u(1,1) * A);
        d(I,1) = B(I,1);
    else
        % Interior points
        S(I,I) = ap;
        b(I,1) = ap;
        for i = 1:n-2
            B(i+1,1) = (rho * u_guess(i,1) * A) - (rho * u_guess(i+1,1) * A);
            d(i+1,1) = B(i+1,1);
            if I == n-1
                % Boundary condition at outlet
                B(I,1) = ((rho * u_guess(i,1) * A) - (rho * u_guess(i+1,1) * A)) + p_cor(end);
                d(I,1) = B(I,1);
            end
        end
    end
end

% Populating upper and lower diagonals for TDMA (Tri-Diagonal Matrix Algorithm)
for I = 1:n-2
    S(I,I+1) = -ae;    % Upper diagonal
    c(I,1) = -ae;
    S(I+1,I) = -aw;    % Lower diagonal
    a(I,1) = -aw;
end

% TDMA Solver for pressure correction
N = length(d);  % Number of rows

% Forward sweep (elimination phase)
c(1) = c(1) / b(1);        % Initial division to avoid zero division
d(1) = d(1) / b(1);        % Normalize first row
for i = 2:N-1
    temp = b(i) - a(i) * c(i-1);
    c(i) = c(i) / temp;     % Update upper diagonal
    d(i) = (d(i) - a(i) * d(i-1)) / temp;  % Update RHS
end

% Backward substitution
d(N) = (d(N) - a(N-1) * d(N-1)) / (b(N) - a(N-1) * c(N-1));
X(N,1) = d(N);
for i = N-1:-1:1
    X(i,1) = d(i,1) - c(i,1) * X(i+1,1);
end

% Updating pressure correction values
for i = 1:n-1
    p_cor(i,1) = X(i,1);
end

% Updating velocities based on pressure correction
for i = 2:n
    for j = 1:n-2
        u(i,1) = u_guess(j,1) + d_co * (p_cor(j) - p_cor(j+1));
    end
end
