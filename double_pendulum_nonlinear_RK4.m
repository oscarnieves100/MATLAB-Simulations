%**************************************************************************
% Program for solving nonlinear double pendulum system using the simple
% RK4 method for a system of 4 1st order equations.
%
% Created by: Oscar A. Nieves
% Last updated: 2017
%**************************************************************************

%% --- Additional Notes
% The following pogram can be used to solve any system of 2 nonlinear (or
% linear) 2nd order differential equations. Simply adjust the initial
% conditions and parameters, and the expressions for the functions fa, fo,
% fX and fv which are the k1 functions of the RK4 algorithm. The program
% also plots the phase portraits for each of the nonlinear systems
% (derivative against its function), so the range of the plots may also
% need to be adjusted to the specific problem. Every other part of the
% program can be left as it stands. 

%% Clear workspace
clc         %clear command prompt
clear       %clear values
close all   %close all figures, tables, diagrams, etc.

%% Define initial conditions, input parameters and constants
m1 = 1;             % mass of first pendulum
m2 = 1.8;           % mass of second pendulum
L1 = 1;             % length of first pendulum
L2 = 1;             % length of second pendulum
g = 9.81;           % gravitational constant
t = 0;              % initial time
X = 1.5;            % initial pangle of second pendulum
a = 0.8;            % initial angle of first pendulum
v = 0;              % initial angular velocity of second pendulum
o = 0;              % initial angular velocity of first pendulum
n = 10000;          % number of points in numerical approximation
h = 50/n;           % step size in the numerical approximation

%% Define variables and initial conditions for animation
X0 = X;
t0 = t;
v0 = v;
a0 = a;
o0 = o;

%% Define variables, vectors, matrices, etc
time = zeros(1,n); 
position = zeros(1,n);
angle = zeros(1,n);
velocity = zeros(1,n);
omega = zeros(1,n);

%% Define functions for RK4 algorithm
fa = @(tf,af,of,Xf,vf) of;
fX = @(tf,af,of,Xf,vf) vf;
fo = @(tf,af,of,Xf,vf) (1-(m2*cos(af-Xf)^2)/(m1+m2))^(-1)*1/(L1*m1+L1*m2)*(m2*L2*(-fa(tf,af,of,Xf,vf)^2*...
    sin(af - Xf) + g/L1*sin(Xf))*cos(af-Xf) - m2*L2*fX(tf,af,of,Xf,vf)*sin(af-Xf) -...
    g*(m1+m2)*sin(af));
fv = @(tf,af,of,Xf,vf) -L1/L2*(fo(tf,af,of,Xf,vf)*cos(af-Xf)-...
    fa(tf,af,of,Xf,vf)^2*sin(af-Xf) + g/L1*sin(Xf));
% --- Note: if there are any coupled equations, uncouple them before typing
% in their expressions.

%% Initiate RK4 method
tic
for i = 1:n
    % place values into vectors for plotting
    time(1,i) = t;
    angle(1,i) = a;
    omega(1,i) = o;
    position(1,i) = X;
    velocity(1,i) = v;
    
    % solve the system for the angular displacement equations
    ka1 = fa(t, a, o, X, v);
    ko1 = fo(t, a, o, X, v);
    kX1 = fX(t, a, o, X, v);
    kv1 = fv(t, a, o, X, v);
    ka2 = fa(t + 0.5*h, a + 0.5*h*ka1, o + 0.5*h*ko1, X + 0.5*h*kX1, v + 0.5*h*kv1);
    ko2 = fo(t + 0.5*h, a + 0.5*h*ka1, o + 0.5*h*ko1, X + 0.5*h*kX1, v + 0.5*h*kv1);
    kX2 = fX(t + 0.5*h, a + 0.5*h*ka1, o + 0.5*h*ko1, X + 0.5*h*kX1, v + 0.5*h*kv1);
    kv2 = fv(t + 0.5*h, a + 0.5*h*ka1, o + 0.5*h*ko1, X + 0.5*h*kX1, v + 0.5*h*kv1);
    ka3 = fa(t + 0.5*h, a + 0.5*h*ka2, o + 0.5*h*ko2, X + 0.5*h*kX2, v + 0.5*h*kv2);
    ko3 = fo(t + 0.5*h, a + 0.5*h*ka2, o + 0.5*h*ko2, X + 0.5*h*kX2, v + 0.5*h*kv2);
    kX3 = fX(t + 0.5*h, a + 0.5*h*ka2, o + 0.5*h*ko2, X + 0.5*h*kX2, v + 0.5*h*kv2);
    kv3 = fv(t + 0.5*h, a + 0.5*h*ka2, o + 0.5*h*ko2, X + 0.5*h*kX2, v + 0.5*h*kv2);
    ka4 = fa(t + h, a + h*ka3, o + h*ko3, X + h*kX3, v + h*kv3);
    ko4 = fo(t + h, a + h*ka3, o + h*ko3, X + h*kX3, v + h*kv3);
    kX4 = fX(t + h, a + h*ka3, o + h*ko3, X + h*kX3, v + h*kv3);
    kv4 = fv(t + h, a + h*ka3, o + h*ko3, X + h*kX3, v + h*kv3);
    
    % Update current values
    a = a + h/6*(ka1 + 2*ka2 + 2*ka3 + ka4);
    o = o + h/6*(ko1 + 2*ko2 + 2*ko3 + ko4);
    X = X + h/6*(kX1 + 2*kX2 + 2*kX3 + kX4);
    v = v + h/6*(kv1 + 2*kv2 + 2*kv3 + kv4);
    t = t + h;
end
toc

%% Plot the resulting functions
plot(time,angle,'r'); xlabel('time (s)'); ylabel('angle (rad)');figure;

plot(time,position,'b'); xlabel('time (s)'); ylabel( ...
    'cart displacement (m)');figure;

plot(time,omega,'r'); xlabel('time (seconds)'); ylabel(...
    'angular velocity (rad/s)');figure;

plot(time,velocity,'b'); xlabel('time (seconds)'); ylabel(...
    'cart velocity (m/s)');

%% Create animated inverted pendulum system
% Define pendulum
pivot = [0 0];                       % pivot position
r = 0.1;                             % radius of pendulum bob
r2 = 0.1;

pos = pivot - (L1*[sin(a0) cos(a0)]);  % Bob 1 position
pos2 = pos - (L2*[sin(X0) cos(X0)]); % Bob 2 position

% --- Rendering parameters
% set scale on x-y plane for animation
figure;
axes = gca;
xlim(axes, [-(L1 + L2 + r + r2) (L1 + L2 + r + r2)]);  
ylim(axes, [-(L1 + L2 + r + r2) (L1 + L2 + r + r2)]);

% create shapes for pendulum and cart
Bob = rectangle('Position', [(pos - r/2) r r], 'Curvature', [1,1], ...
    'FaceColor', 'r');                                        % bob render
hold on
rod = line([pivot(1) pos(1)], [pivot(2) pos(2)]);             % rod render
hold on
Bob2 = rectangle('Position', [(pos2 - r2/2) r2 r2], ...
     'Curvature', [1,1], 'FaceColor', 'b');    % cart render
hold on
rod2 = line([pos(1) pos2(1)], [pos(2) pos2(2)]);
hold off

% --- Initiate animation
for time1 = 1:n                % set same time scale as for solution to ODEs
    drawnow;                  % draw elements
    X0 = position(1,time1);
    a0 = angle(1,time1);
    
    % update position of bob, rod and cart
    pos = pivot - (L1*[sin(a0) cos(a0)]);
    pos2 = pos - (L2*[sin(X0) cos(X0)]);
    set(Bob, 'Position',[(pos - r/2) r r]);
    set(Bob2, 'Position', [(pos2 - r2/2) r2 r2]);
    set(rod, 'XData', [pivot(1) pos(1)], 'YData', [pivot(2) pos(2)]);
    set(rod2, 'XData', [pos(1) pos2(1)], 'YData', [pos(2) pos2(2)]);
end 