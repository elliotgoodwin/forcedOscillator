%%%%%%%%%%%%%%%%%%%%%%%%%%    Project 2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elliot Goodwin - 9621958
% Simulates SHM using four different numerical methods.  The best method,
% (compared to the analytical solution when F(t) = 0) will be used to
% solve F(t) =/ 0 (no analytical solution)


clear all;
close all;
format long;


%%%%%%%%%%%%%%%%% Input Parameters and Initial Conditions %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initially set the force term to 0
params2 = 0;

% open user interface for variable input
prompt = {'Spring constant  [N/m]:', 'Mass  [kg]:',...
    'Damping coefficient  [kg/s]:', 'Initial position  [m]:',...
    'Initial velocity  [m/s]:', 'Simulation time  [s]:', 'Step size  [s]'};
defAns = {'1.73', '3.03', '0.1', '0', '-1', '100', '0.1'};
params = inputdlg(prompt, 'Input Parameters', 1, defAns);

% if user presses cancel, exit programme
if isempty(params) == 1
    clear all;
    error('User terminated programme.');
end

k = str2double(params{1, 1});
m = str2double(params{2, 1});
b = str2double(params{3, 1});
x0 = str2double(params{4, 1});
v0 = str2double(params{5, 1});
T = str2double(params{6, 1});
h = str2double(params{7, 1});

% check that input parameters are sensible
% loop whilst parameters are invalid
while k <= 0 | b < 0 | m <= 0 | T <= 0 | h <= 0;

    % display help window
    uiwait(helpdlg({'Please ensure parameters have appropriate values:',...
        '', 'Spring constant > 0', 'Mass > 0', 'Damping coefficient => 0',...
        'Simulation run time > 0', 'Time step > 0'}, 'Error'));
    % re-ask for parameter input
    params = inputdlg(prompt, 'Input Parameters', 1, defAns);
   
    k = str2double(params{1, 1});
    m = str2double(params{2, 1});
    b = str2double(params{3, 1});
    x0 = str2double(params{4, 1});
    v0 = str2double(params{5, 1});
    T = str2double(params{6, 1});
    h = str2double(params{7, 1});
    
    % is user presses cancel, exit programme
    if isempty(params) == 1
        error('User terminated programme.');
    end

end
    

% create an array of t values
% 0 <= t <= T, in increments of h
t = [0:h:(T-h)];
t = t';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Analytical Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = k/m;
gamma = b/m;
arg = w0 - (gamma/2)^2;
w = sqrt(arg);

prefactor = exp( - ((gamma.*t)/2) );
C_1 = ( v0 + (gamma.*x0)./2 ) ./ w;
C_2 = x0;

% light damping
if arg > 0
    x = prefactor.*( C_1.*sin(w.*t) + C_2.*cos(w.*t) );
    v = -(gamma./2).*prefactor.*( C_1.*sin(w.*t) + C_2.*cos(w.*t) )...
        + prefactor.*( w.*C_1.*cos(w.*t) - w.*C_2.*sin(w.*t) );
    E = 0.5.*k.*x.^2 + 0.5.*m.*v.^2;
% heavy damping
elseif arg < 0 
    x = prefactor.*( C_1.*sinh(w.*t) + C_2.*cosh(w.*t) );
    v = -(gamma./2).*prefactor.*( C_1.*sinh(w.*t) + C_2.*cosh(w.*t) )...
        + prefactor.*( w.*C_1.*cosh(w.*t) + w.*C_2.*sinh(w.*t) );
    E = 0.5.*k.*x.^2 + 0.5.*m.*v.^2;
% critical damping
else 
    x = C_1.*w.*t + C_2;
    v = C_1.*w;
    E = 0.5.*k.*x.^2 + 0.5.*m.*v.^2;
end

% plot each method on the same axis
figure(1);
% Euler's Method
[x_euler, v_euler, E_euler] = eulerfunc(x0, v0, T, k, b, m, h, params2);
plot(t, x_euler, 'y-');     hold on;
% Improved Euler's Method
[x_improvedeuler, v_improvedeuler, E_improvedeuler] = improvedeuler(x0,...
    v0, T, k, b, m, h, params2);
plot(t, x_improvedeuler, 'r-');     hold on;
% Verlet Method
[x_verlet, v_verlet, E_verlet] = verletfunc(x0, v0, T, k, b, m, h, params2);
plot(t, x_verlet, 'b-');     hold on;
% Euler-Cromer Method
[x_cromer, v_cromer, E_cromer] = eulercromer(x0, v0, T, k, b, m, h, params2);
plot(t, x_cromer, 'g');   hold on;
% Analytical solution
plot(t, x, 'k-'); hold on;

%title('x(t) Evaluated Using Four Different Integrators');   hold on;
legend('Euler''s Method', 'Improved Euler Method', 'Verlet Method',...
    'Euler-Cromer Method', 'Analytical Solution', 'Location',...
    'Southwest');    hold off;
xlabel('Time  [s]');
ylabel('Amplitude of oscillation  [m]');
xlim([0, 70]);

% find 'residuals' - difference between numerical and analytic solutions
% for each method
y_euler = x_euler - real(x);
y_improvedeuler = x_improvedeuler - real(x);
y_verlet = x_verlet - real(x);
y_cromer = x_cromer - real(x);
% plot 'residuals'
figure(2);
plot(t, y_euler, 'y');   hold on;
plot(t, y_improvedeuler, 'r'); hold on;
plot(t, y_verlet, 'b');  hold on;
plot(t, y_cromer, 'g');  hold on;
% this plot is just for reference
plot(t, zeros(size(x)), 'k'); hold on;
%title('Residuals'); hold on;
legend('Euler''s Method', 'Improved Euler Method', 'Verlet Method',...
    'Euler-Cromer Method', 'Analytical Solution', 'Location',...
    'best');    hold off;
xlabel('Time  [s]');
ylabel('x_{numeric} - x_{analytic}  [m]');
ylim([-0.1, 0.1]);

% find normal of residuals to quantify how good each method
normr_euler = sum((y_euler.^2));
normr_improvedeuler = sum(y_improvedeuler.^2);
normr_verlet = sum(y_verlet.^2);
normr_cromer = sum(y_cromer.^2);
snormr = [normr_euler, normr_improvedeuler, normr_verlet, normr_cromer];



%%%%%%%%%%%%%%%%%%%%%%% Write Data to File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% write euler data out to file
% create a structure containing information to be written
A = [t, x_euler, v_euler, E_euler];
% open a file for data to be written to
eulerdata = fopen('euler_data.txt','w');
% check file is open
if eulerdata < 0 
    error('Cannot write to file.');
else
    % if file is open, write data to .txt file
    fprintf(eulerdata, 'Euler''s method:\r\n\r\n');
    fprintf(eulerdata,'%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]',...
        'v  [m/s]', 'E  [J]');
    fprintf(eulerdata,'%4.2f %12.4f %12.4f %12.4f \r\n', A');
    fclose(eulerdata);
end


% write improved euler data out to file
% create a structure containing information to be written
B = [t, x_improvedeuler, v_improvedeuler, E_improvedeuler];
% open a .txt file for data to be written to
improvedeulerdata = fopen('improvedeuler_data.txt','w');
% check file is open
if improvedeulerdata < 0 
    error('Cannot write to file.');
else
    % if file is open, write data to .txt file
    fprintf(improvedeulerdata,...
        'Improved Euler''s method:\r\n\r\n');
    fprintf(improvedeulerdata,'%5s %10s %14s %10s\r\n', 't  [s]',...
        'x  [m]', 'v  [m/s]', 'E  [J]');
    fprintf(improvedeulerdata,'%4.2f %12.4f %12.4f %12.4f \r\n', B');
    fclose(improvedeulerdata);
end


% add verlet data to a structure to be written to .txt file
C = [t, x_verlet, v_verlet, E_verlet];
% open .txt file
verletdata = fopen('verlet_data.txt','w');
% check file is open
if verletdata < 0 
    error('Cannot write to file.');
else
    % if file is open, write data to .txt file
    fprintf(verletdata, 'Verlet method:\r\n\r\n');
    fprintf(verletdata,'%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]',...
        'v  [m/s]', 'E  [J]');
    fprintf(verletdata,'%4.2f %12.4f %12.4f %12.4f \r\n', C');
    fclose(verletdata);
end


% add euler-cromer data to structure to be written to .txt file
D = [t, x_cromer, v_cromer, E_cromer]; 
% open .txt file
cromerdata = fopen('cromer_data.txt','w');
% check file is open
if cromerdata < 0 
   error('Cannot write to file.');
else
    % if file is open, write data to .txt file
    fprintf(cromerdata, 'Euler-Cromer method:\r\n\r\n');
    fprintf(cromerdata, '%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]',...
        'v  [m/s]', 'E  [J]');
    fprintf(cromerdata,'%4.2f %12.4f %12.4f %12.4f \r\n', D');
    fclose(cromerdata);
end



%%%%%%%%%%%%%%%%%%%% Investigating Time Step Size %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% use Euler's method with two different values of step size
% b = 0 (no damping) ==> analytical solution has constant energy
% create time arrays for each step size
t_small = [0:0.001:(T-0.001)];
t_med = [0:0.01:(T-0.01)];
t_large = [0:0.1:(T-0.1)];

figure(3);
[x_small, v_small, E_small] = eulerfunc(x0, v0, T, k, 0, m, 0.001, params2);
plot(t_small, x_small, 'g-');     hold on;
[x_med, v_med, E_med] = eulerfunc(x0, v0, T, k, 0, m, 0.01, params2);
plot(t_med, x_med, 'r-');     hold on;
[x_large, v_large, E_large] = eulerfunc(x0, v0, T, k, 0, m, 0.1, params2);
plot(t_large, x_large, 'b-');     hold on;
plot(t, x, 'k-');   hold off;
%title('Euler''s Method with varying step size');
legend('Step size = 0.001 s', 'Step size = 0.01 s', 'Step size = 0.1 s',...
    'Analytical solution', 'Location', 'best');
xlabel('Time  [s]');
ylabel('Amplitude of oscillation  [m]');
xlim([0, 50]);

% plot E(t) for each step size to compare to analytical solution
figure(4);
plot(t_small, E_small, 'g');     hold on;
plot(t_med, E_med, 'r');         hold on;
plot(t_large, E_large, 'b');     hold on;
%title('E(t) calculated using Euler''s Method for varying step size');
legend('Step size = 0.001 s', 'Step size = 0.01 s', 'Step size = 0.1 s',...
    'Location', 'best');
xlabel('Time  [s]');
ylabel('Energy  [J]');
ylim([0, 10]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read in x and v data from .txt file
temp = importdata('euler_data.txt');
data = temp.data;
t_euler = data(:, 1);
x_euler = data(:, 2);
v_euler = data(:, 3);
E_euler = data(:, 4);
clear temp data;

temp = importdata('improvedeuler_data.txt');
data = temp.data;
t_improvedeuler = data(:, 1);
x_improvedeuler = data(:, 2);
v_improvedeuler = data(:, 3);
E_improvedeuler = data(:, 4);
clear temp data;

temp = importdata('verlet_data.txt');
data = temp.data;
t_verlet = data(:, 1);
x_verlet = data(:, 2);
v_verlet = data(:, 3);
E_verlet = data(:, 4);
clear temp data;

temp = importdata('cromer_data.txt');
data = temp.data;
t_cromer = data(:, 1);
x_cromer = data(:, 2);
v_cromer = data(:, 3);
E_cromer = data(:, 4);
clear temp data;


% plot E(t) to compare methods
% E(t) should be have gradient = 0 when b = 0
figure(5);
%title('E(t) for four different integrators');   hold on;
plot(t, E_euler, 'y');   hold on;
plot(t, E_improvedeuler, 'r');   hold on;
plot(t, E_verlet, 'b');  hold on;
plot(t, E_cromer, 'g');  hold on;
plot(t, E, 'k'); hold off;
legend('Euler''s Method', 'Improved Euler Method', 'Verlet Method',...
    'Euler-Cromer Method', 'Analytic Solution', 'Location', 'best');    hold off;
xlabel('Time  [s]');
ylabel('Energy [J]');
ylim([0, 5]);



%%%%%%%%%%%%%%%%%%%%%%%%%%% Phase Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% phase space plots, x(v)
figure(6);
plot(x_euler, v_euler); hold on;
plot(x, v); hold off;
%title('x(v) Euler method');
xlabel('Amplitude of oscillation  [m]');
ylabel('Instantaneous velocity  [ms^-^1]');
legend('Euler''s method', 'Analytical solution', 'Location', 'best');

figure(7);
plot(x_improvedeuler, v_improvedeuler); hold on;
plot(x, v); hold off;
%title('x(v) Improved Euler method');
xlabel('Amplitude of oscillation  [m]');
ylabel('Instantaneous velocity  [ms^-^1]');
legend('Improved Euler''s method', 'Analytical solution', 'Location',...
    'best');

figure(8);
plot(x_verlet, v_verlet);   hold on;
plot(x, v); hold off;
%title('x(v) Verlet method');
xlabel('Amplitude of oscillation  [m]');
ylabel('Instantaneous velocity  [ms^-^1]');
legend('Verlet method', 'Analytical solution', 'Location', 'best');

figure(9);
plot(x_cromer, v_cromer);   hold on;
plot(x, v); hold off;
%title('x(v) Euler-Cromer method');
xlabel('Amplitude of oscillation  [m]');
ylabel('Instantaneous velocity  [ms^-^1]');
legend('Euler-Cromer method', 'Analytical solution', 'Location', 'best');



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Critical Damping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% All further analysis performed using Verlet Method
% Investigate critical damping
b_crit = 2*sqrt(k*m);
figure(10);
[x_crit, v_crit, E_crit] = verletfunc(x0, v0, T, k, b_crit, m, h, params2);
plot(t, x_crit, 'b-');     hold on;
[x_2crit, v_2crit, E_2crit] = verletfunc(x0, v0, T, k, 2*b_crit,...
    m, h, params2);
plot(t, x_2crit, 'r-');     hold on;
[x_halfcrit, v_halfcrit, E_halfcrit] = verletfunc(x0, v0, T, k,...
    0.5*b_crit, m, h, params2);
plot(t, x_halfcrit, 'g-');     hold off;
%title('Investigating critical damping');
legend('b = 2\surd(km)', 'b = 4\surd(km)', 'b = \surd(km)', 'Location',...
    'best');
xlabel('Time  [s]');
ylabel('Amplitude of oscillations  [m]');
xlim([0, 30]);



%%%%%%%%%%%%%%%%%%%%%%%%% Investingating Impulse %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create a few force vectors, with different magnitudes, directions and
% times of application
F1 = ones(size(t));
F1(100) = 10;
F2 = ones(size(t));
F2(650) = -4;
F3 = ones(size(t));
F3(370) = -6;

figure(11);
[x_F1, v_F1, E_F1] = verletfunc(x0, v0, T, k, b, m, h, F1);
plot(t, x_F1, 'b-');     hold on;
[x_F2, v_F2, E_F2] = verletfunc(x0, v0, T, k, b, m, h, F2);
plot(t, x_F2, 'r-');     hold on;
[x_F3, v_F3, E_F3] = verletfunc(x0, v0, T, k, b, m, h, F3);
plot(t, x_F3, 'g-');     hold off;
%title('x(t) for an instantaneously forced oscillator');
legend('F = 10 N, t = 10 s', 'F = -4 N, t = 65 s', 'F = -6 N, t = 37 s',...
    'Location', 'best');
xlabel('Time  [s]');
ylabel('Amplitude of oscillation  [m]');


% create sinusoidal force vectors
Fsine2 = ones(size(t)).*sin(2.*w.*t);
Fsine3 = ones(size(t)).*sin(0.5.*w.*t);

figure(12);
[x_Fsine2, v_Fsine2, E_Fsine2] = verletfunc(x0, v0, T, k, 0, m, h, Fsine2);
plot(t, x_Fsine2, 'r-');     hold on;
[x_Fsine3, v_Fsine3, E_Fsine3] = verletfunc(x0, v0, T, k, 0, m, h, Fsine3);
plot(t, x_Fsine3, 'b-');     hold on;
%title('x(t) for a sinusoidally forced oscillator');
xlabel('Time  [s]');
ylabel('Amplitude of oscillation  [m]');
legend('F = sin(2\omegat)', 'F = sin(0.5\omegat)', 'Location', 'best');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Resonance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for n = 1:200;
    
    u(n) = n.*w./100;
    F = sin(u(n)*t);
    x1 = abs(verletfunc(x0, v0, T, k, 0, m, h, F));
    x_amp1(n) = max(x1);
    x2 = abs(verletfunc(x0, v0, T, k, 0.01, m, h, F));
    x_amp2(n) = max(x2);
    x3 = abs(verletfunc(x0, v0, T, k, 0.1, m, h, F));
    x_amp3(n) = max(x3);
    x4 = abs(verletfunc(x0, v0, T, k, 0.5, m, h, F));
    x_amp4(n) = max(x4);
    x5 = abs(verletfunc(x0, v0, T, k, 1, m, h, F));
    x_amp5(n) = max(x5);
    x6 = abs(verletfunc(x0, v0, T, k, b_crit, m, h, F));
    x_amp6(n) = max(x6);
    
    
end

figure(13);
plot(u, x_amp1);    hold on;
plot(u, x_amp2);    hold on;
plot(u, x_amp3);    hold on;
plot(u, x_amp4);    hold on;
plot(u, x_amp5);    hold on;
plot(u, x_amp6);    hold off;
%title('Resonance curve');
xlabel('\omega''  [rads^-^1]');
ylabel('Maximum amplitude of oscillation  [m]');
legend('b = 0 kgs^-^1', 'b = 0.01 kgs^-^1', 'b = 0.1 kgs^-^1',...
     'b = 0.5 kgs^-^1', 'b = 1 kgs^-^1', 'b = 2\surd(km)', 'Location', 'best')
xlim([0, 1.6]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


saveas(figure(1),'x(t) four integrators','png');
saveas(figure(2),'residuals','png');
saveas(figure(3),'x(t) euler step size','png');
saveas(figure(4),'E(t) euler step size','png');
saveas(figure(5),'E(t) four integrators','png');
saveas(figure(6),'x(v) euler','png');
saveas(figure(7),'x(v) improved euler','png');
saveas(figure(8),'x(v) verlet','png');
saveas(figure(9),'x(v) euler-cromer','png');
saveas(figure(10),'critical damping','png');
saveas(figure(11),'instantaneous forcing','png');
saveas(figure(12),'sinusoidal forcing','png');
saveas(figure(13),'resonance curve','png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

