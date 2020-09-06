%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares-Analysis (LSA) for linear fit
%
% Description: generates a random data-set around a linear curve with
% preset gradient and then uses LSA matrix methods to find the coefficients
% of the line of best-fit, and overlays the fit-line to the dataset:
% Equation --> y = b0 + b1*x
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Generate dataset
x = 0:0.1:10;
S = 2*x + randn(size(x));

% Use LSA to find line of best fit
n = length(S);
sum_x = sum(x);
sum_x2 = sum(x.^2);
sum_S = sum(S);
sum_xS = sum(x.*S);

A = [ [n, sum_x]; ...
     [sum_x, sum_x2] ];
RHS = [sum_S, sum_xS].';
b = A\RHS;
b0 = b(1);
b1 = b(2);

% Fit straight line to data and generate plots
y = b0 + b1*x;

figure(1);
set(gcf,'color','w');
scatter(x,S); hold on;
plot(x,y,'r','LineWidth',3); hold off;
legend('Raw Data','Fit-line'); legend boxoff;
legend('Location','northwest');
title(['b_0 = ' num2str(round(b0,2)) ', b_1 = ' num2str(round(b1,2))]);
xlabel('x'); 
ylabel('y');
axis tight;
set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%