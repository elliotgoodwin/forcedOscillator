function [ x, v, E ] = improveuler(x0, v0, T, k, b, m, h, F)
% Improved Euler Method
%   Improves Euler's method, by adding an extra term from the taylor
%   expansion

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
a = zeros(z, 1);
F = zeros(z, 1);
E = zeros(z, 1);

% initialise variables from function input arguments
x(1) = x0;
v(1) = v0;
a(1) = -(b./m).*v(1) - (k./m)*x(1) + F(1)./m;
E(1) = 0.5*m*v(1).^2 + 0.5*k*x(1).^2;

% find next value of x
x(2) = x(1) + h.*v(1) + ((h.^2)./2).*a(1);

% loop z times, to find z values of x
for n = 2:z;
   
    v(n) = v(n-1) + h.*a(n-1);
    x(n) = x(n-1) + h.*v(n-1) + ((h.^2)./2).*a(n-1);
    a(n) = -(b./m).*v(n) - (k./m).*x(n);
    E(n) = 0.5*m*v(n).^2 + 0.5*k*x(n).^2;
    
end;

end

