%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo integration (naive MC) in 2 dimensions. It sets up a volume
% around the domain defined by the user, and then places uniformly
% generated random numbers at coordinate (x,y,z) (each coordinate sampled
% separately). If the random point p(x0,y0,z0) falls within the volume of
% integration between the bounding surfaces, then it counts it as a 'hit'
% and adds one to the counter, otherwise it ignores it. The integral is
% approximated by taking the ratio of hits to total points generated.
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Inputs
iterations = 1000;
N = 1000;
xlow = -1;
xhigh = 1;
ylow = @(x) x.^2 - 1;
yhigh = @(x) 1 - x.^2;
z = @(x,y) exp(-x.^2 - y.^2); %surface

% Volume of rectangular prism
xv = linspace(xlow,xhigh,N);
xmax = xhigh; xmin = xlow;
ylmin = unique(min(ylow(xv))); ylmax = unique(max(ylow(xv)));
yhmin = unique(min(yhigh(xv))); yhmax = unique(max(yhigh(xv)));
ymin = min([ylmin yhmin]); ymax = max([ylmax yhmax]);
yv = linspace(ymin, ymax, N);
[X,Y] = meshgrid(xv,yv);
ZV = z(X,Y);
zmax = unique(max(max(ZV))); zmin = unique(min(min(ZV)));
if (zmin < 0) && (zmax > 0)
    rect = abs(zmax - zmin)*(xmax - xmin)*(ymax - ymin);
elseif (zmax <= 0)
    rect = abs(zmin)*(xmax - xmin)*(ymax - ymin);
elseif (zmin >= 0)
    rect = zmax*(xmax - xmin)*(ymax - ymin);
end

% Monte Carlo simulation
random_x = rand(iterations,N);
random_y = rand(iterations,N);
random_z = rand(iterations,N);
xrand = xmin + (xmax - xmin)*random_x;
yrand = ymin + (ymax - ymin)*random_y;
yflow = ylow(xrand);
yfhigh = yhigh(xrand);
yfL = min(yflow,yfhigh);
yfH = max(yflow,yfhigh);
if (zmin < 0) && (zmax > 0)
    zrand = zmin + (zmax - zmin)*random_z;
elseif (zmax <= 0)
    zrand = zmin*random_z;
elseif (zmin >= 0)
    zrand = zmax*random_z;
end
zf = z(xrand,yrand);
H_pos = sum(((yrand > yfL) & (yrand < yfH)) & ...
    ((zf > 0) & (zrand > 0) & (zrand < zf)),2);
H_neg = sum(((yrand > yfL) & (yrand < yfH)) & ...
    ((zf < 0) & (zrand < 0) & (zrand > zf)),2);
I = mean((H_pos - H_neg)/N*rect);

% Exact value (when x = [xlow, xhigh])
exact_I = integral(@(x) 0.5*exp(-x.^2).*sqrt(pi).*...
    (erf(1-x) + erf(1-x.^2)), xlow, xhigh);
error_I = 100*abs( (exact_I - I)/exact_I );

% Display results
disp(['Monte Carlo integral = ' num2str(I)]);
disp(['Exact integral = ' num2str(exact_I)]);
disp(['%error = ' num2str(error_I)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%