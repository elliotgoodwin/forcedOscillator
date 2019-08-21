function [ x, v, E ] = verlet(x0, v0, T, k, b, m, h, f)
% Verlet Method
%   Uses the Verlet method to approximate x(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input arguments:  - x0 = initial position
%                   - v0 = initial velocity
%                   - T = total simulation time
%                   - k = spring constant
%                   - b = damping coefficient
%                   - m = mass
%                   - h = time step increment
%                   - F = Force term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output arguments  - t = time
%                   - x = position
%                   - E = total energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% number of data points, z
z = T./h;
% initalise x, v arrays with z cells
x = zeros(z, 1);
v = zeros(z, 1);
F = ones(z, 1).*f;
E = zeros(z, 1);

% initialise variables from function input arguments
x(1) = x0;
v(1) = v0;

E(1) = 0.5*m*v(1).^2 + 0.5*k*x(1).^2;

D = (2*m + b*h);
Anum = 2*(2*m - k*(h^2)); 
Bnum = (b*h - 2*m);

A = Anum/D;
B = Bnum/D;
C = 2*h/D;

% use improved euler method to find next value of x
x(2) = x(1)*(1 - ((h^2)*k)/(2*m)) + h*v(1)*(1 - ((h^2)*b)/(2*m));

% loop verlet method z times
for n = 2:z;
    
    x(n+1) = A*x(n) + B*x(n-1) + C*F(n);
    % find verlet velocity
    v(n) = (x(n+1) - x(n-1))./(2*h);
    % calculate energy at each point in time
    E(n) = 0.5*m*v(n).^2 + 0.5*k*x(n).^2;
    
end

% remove last cell of array to be consistent with other arrays
x = x(1:end-1);

end

