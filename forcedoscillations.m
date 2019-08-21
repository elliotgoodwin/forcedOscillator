%%%%%%%%%%%%%%%%%%%%%%%%%%  forcedOscillations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates SHM using four different numerical methods.  The best method,
% (compared to the analytical solution when F(t) = 0) will be used to
% solve F(t) =/= 0 (no analytical solution)

clear all;
close all;
format long;

% open user interface for variable input
prompt = {'Spring constant  [N/m]:', 'Mass  [kg]:',...
    'Damping coefficient  [kg/s]:', 'Initial position  [m]:',...
    'Initial velocity  [m/s]:', 'Simulation time  [s]:', 'Step size  [s]'};
defAns = {'1.73', '3.03', '0.1', '1', '0', '100', '0.1'};
params = inputdlg(prompt, 'Input Parameters', 1, defAns);

% if user presses cancel, exit programme
if isempty(params) == 1
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
err = 'Please ensure all parameters are real and positive.';
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

% plot each method on the same axis
figure(1);
title('x(t) evaluated using four integrators');
% Euler's Method
[x_euler, v_euler, E_euler] = eulerfunc(x0, v0, T, k, b, m, h);
plot(t, x_euler, 'g-');     hold on;
% Improved Euler's Method
[x_improvedeuler, v_improvedeuler, E_improvedeuler] = improvedeuler(x0, v0, T, k, b, m, h);
plot(t, x_improvedeuler, 'r-');     hold on;
% Verlet Method
[x_verlet, v_verlet, E_verlet] = verletfunc(x0, v0, T, k, b, m, h);
plot(t, x_verlet, 'b-');     hold on;
% Euler-Cromer Method
[x_cromer, v_cromer, E_cromer] = eulercromer(x0, v0, T, k, b, m, h);
plot(t, x_cromer);   hold on;


% Analytical solution
% nb analytical solution is wrong, probably the amplitude bit
omega = sqrt(k./m);
gamma = b./m;
arg = (i.*omega - gamma./2).*t;
x = x0.*exp(arg);
plot(t, x, 'k-'); hold off;


% find 'residuals' - difference between numerical and analytic solutions
% for each method
y_euler = x_euler - real(x);
y_improvedeuler = x_improvedeuler - real(x);
y_verlet = x_verlet - real(x);
y_cromer = x_cromer - real(x);
% plot 'residuals'
figure(2);
title('Residuals');
plot(t, y_euler, 'g');   hold on;
plot(t, y_improvedeuler, 'r'); hold on;
plot(t, y_verlet, 'b');  hold on;
plot(t, y_cromer);  hold on;
% this plot is just for reference
plot(t, zeros(size(x))); hold off;

% find normal of residuals to quantify how good each method
normr_euler = sum((y_euler.^2));
normr_improvedeuler = sum(y_improvedeuler.^2);
normr_verlet = sum(y_verlet.^2);
normr_cromer = sum(y_cromer.^2);
snormr = [normr_euler, normr_improvedeuler, normr_verlet, normr_cromer];



%%%%%%%%%%%%%%%%%%%%%%% write data to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    fprintf(eulerdata,'%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]', 'v  [m/s]', 'E  [J]');
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
    fprintf(improvedeulerdata, 'Improved Euler''s method to 12 dp precision:\r\n\r\n');
    fprintf(improvedeulerdata,'%5s %10s %22s %18s\r\n', 't  [s]', 'x  [m]', 'v  [m/s]', 'E  [J]');
    fprintf(improvedeulerdata,'%4.2f %20.12f %20.12f %20.12f \r\n', B');
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
    fprintf(verletdata,'%5s %8s %20s %16s\r\n', 't  [s]', 'x  [m]', 'v  [m/s]', 'E  [J]');
    fprintf(verletdata,'%4.2f %18.12f %18.12f %18.12f \r\n', C');
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
    fprintf(cromerdata,'%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]', 'v  [m/s]', 'E  [J]');
    fprintf(cromerdata,'%4.2f %12.4f %12.4f %12.4f \r\n', D');
    fclose(cromerdata);
end



%%%%%%%%%%%%%%%%%%%%%% example of rounding errors %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the following file is used to illustrate the effect of rounding errors in
% calculations

% repeat write file for improved euler data with only 4 dp precision
E = [t, x_improvedeuler, v_improvedeuler, E_improvedeuler];
roundingexample = fopen('improvedeuler_data_4dp.txt','w');
errmsg = 'Error: cannot write to file.';
if roundingexample < 0 
    error('Cannot write to file.');
else
    fprintf(roundingexample, 'Improved Euler''s method to 4 dp precision:\r\n\r\n');
    fprintf(roundingexample, '%5s %10s %14s %10s\r\n', 't  [s]', 'x  [m]', 'v  [m/s]', 'E  [J]');
    fprintf(roundingexample, '%4.2f %12.4f %12.4f %12.4f \r\n', B');
    fclose(roundingexample);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% read in x and v data from .txt file
temp = importdata('euler_data.txt');
data = temp.data;
t_euler1 = data(:, 1);
x_euler1 = data(:, 2);
v_euler1 = data(:, 3);
E_euler1 = data(:, 4);
clear temp data;

temp = importdata('improvedeuler_data.txt');
data = temp.data;
t_improvedeuler1 = data(:, 1);
x_improvedeuler1 = data(:, 2);
v_improvedeuler1 = data(:, 3);
E_improvedeuler1 = data(:, 4);
clear temp data;

temp = importdata('verlet_data.txt');
data = temp.data;
t_verlet1 = data(:, 1);
x_verlet1 = data(:, 2);
v_verlet1 = data(:, 3);
E_verlet1 = data(:, 4);
clear temp data;

temp = importdata('cromer_data.txt');
data = temp.data;
t_cromer1 = data(:, 1);
x_cromer1 = data(:, 2);
v_cromer1 = data(:, 3);
E_cromer1 = data(:, 4);
clear temp data;


% plot E(t) to compare methods
% E(t) should be have gradient = 0 when b = 0
figure(3);
title('E(t) for four different integrators');
plot(t, E_euler1, 'g');   hold on;
plot(t, E_improvedeuler1, 'r');   hold on;
plot(t, E_verlet1, 'b');  hold on;
plot(t, E_cromer1);  hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%% phase plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% phase space plots
figure(4);
title('x(v) euler');
plot(x_euler, v_euler); hold on;
plot(x_euler1, v_euler1);   hold off;
figure(5);
title('x(v) improved euler');
plot(x_improvedeuler, v_improvedeuler); hold on;
plot(x_improvedeuler1, v_improvedeuler1);  hold off;
figure(6);
title('x(v) verlet');
plot(x_verlet, v_verlet);   hold on;
plot(x_verlet1, v_verlet1); hold off;
figure(7);
title('x(v) euler-cromer');
plot(x_cromer, v_cromer, x_cromer1, v_cromer1);


% formula to evaluate differences in data
%diff = B(:, 2) - x_improvedeuler;
%disp(max(diff));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% verlet data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% investigate critical damping
b_crit = 2*sqrt(k*m);
figure(8);
title('investigating critical damping');
[x_crit, v_crit, E_crit] = verletfunc(x0, v0, T, k, b_crit, m, h);
plot(t, x_crit, 'b-');     hold on;
[x_2crit, v_2crit, E_2crit] = verletfunc(x0, v0, T, k, 2*b_crit, m, h);
plot(t, x_2crit, 'r-');     hold on;
[x_halfcrit, v_halfcrit, E_halfcrit] = verletfunc(x0, v0, T, k, 0.5*b_crit, m, h);
plot(t, x_halfcrit, 'g-');     hold on;

disp('done');
