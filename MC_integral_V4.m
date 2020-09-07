%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo integration (naive MC) in 1 dimension
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Inputs
iterations = 1000;
N = 1000;
a = -pi;
b = pi;
y = @(x) exp(x); % test function

% Area of rectangle
xv = linspace(a,b,N);
yv = y(xv);
ymax = unique(max(yv));
ymin = unique(min(yv));
if (ymin < 0) && (ymax > 0)
    rect = abs(ymax - ymin)*(b - a);
elseif (ymax <= 0)
    rect = abs(ymin)*(b - a);
elseif (ymin >= 0)
    rect = ymax*(b - a);
end

% Monte Carlo simulation
random_x = rand(iterations,N);
random_y = rand(iterations,N);
xrand = a + (b - a)*random_x;
if (ymin < 0) && (ymax > 0)
    yrand = ymin + (ymax - ymin)*random_y;
elseif (ymax <= 0)
    yrand = ymin*random_y;
elseif (ymin >= 0)
    yrand = ymax*random_y;
end
yf = y(xrand);
H_pos = sum((yf > 0) & (yrand > 0) & (yrand < yf),2);
H_neg = sum((yf < 0) & (yrand < 0) & (yrand > yf),2);
I = mean((H_pos - H_neg)/N*rect);

% Exact integral value
exact_I = integral(@(x) y(x), a, b);
error_I = 100*abs( (exact_I - I)/exact_I );

% Display results
disp(['Monte Carlo integral = ' num2str(I)]);
disp(['Exact integral = ' num2str(exact_I)]);
disp(['%error = ' num2str(error_I)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%