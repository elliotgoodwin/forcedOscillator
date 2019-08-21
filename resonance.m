% resonance
close all;
clear all;
x0 = 0;
v0 = -1;
T = 100;
k = 1.73;
b = 0;
m = 3.03;
h = 0.1;
w = k/m;

t = [0:h:(T-h)];
t = t';

for n = 1:40;
    
    u(n) = w - n.*0.01;
    F = sin(u(n)*t);
    x = verletfunc(x0, v0, T, k, 0, m, h, F);
    x_amp(n) = max(x);
    
end

figure(13);
plot(u, x_amp);
title('Resonance curve');
xlabel('\omega  [rads^-^1]');
ylabel('Amplitude of oscillations  [m]');