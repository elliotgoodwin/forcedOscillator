function [ x, v, E ] = eulerfunc(x0, v0, T, k, b, m, h, F)
% Euler Method
%   Uses Euler's method to find x(t), given a number of initial conditions

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
z = T/h;
% initalise x, v, E arrays with z cells
x = zeros(z, 1);
v = zeros(z, 1);
a = zeros(z, 1);
F = zeros(z, 1);
E = zeros(z, 1);
% supply a force of 100 N a t = T/2h
%F(z/2) = -100;

% input initial conditions to x, v arrays
x(1) = x0;
v(1) = v0;
a(1) = -(b./m).*v(1) - (k./m)*x(1) + F(1)./m;
E(1) = 0.5*m*v(1).^2 + 0.5*k*x(1).^2;

% find next value of x
x(2) = x(1) + h*v(1);

% loop z times, to find z values of x
for n = 2:z;
   
    
    v(n) = v(n-1) + h.*a(n-1);
    x(n) = x(n-1) + h*v(n-1);
    a(n) = -(b./m).*v(n) - (k./m)*x(n) + F(n)./m;
    E(n) = 0.5*m*v(n-1).^2 + 0.5*k*x(n-1).^2;
    
end

end
