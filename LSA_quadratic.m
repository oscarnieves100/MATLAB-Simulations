%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares-Analysis (LSA) for linear fit
%
% Description: generates a random data-set around a quadratic curve with
% preset values and then uses LSA matrix methods to find the coefficients
% of the line of best-fit, and overlays the fit-line to the dataset:
% Equation --> y = b0 + b1*x + b2*x^2
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Generate dataset
x = 0:0.1:5;
S = 3*x.^2 - x + 5*randn(size(x));

%% Quadratic LSA
n = length(S);
sum_x = sum(x);
sum_x2 = sum(x.^2);
sum_x3 = sum(x.^3);
sum_x4 = sum(x.^4);
sum_S = sum(S);
sum_xS = sum(x.*S);
sum_x2S = sum(x.^2.*S);

A = [ [n, sum_x, sum_x2]; ...
      [sum_x, sum_x2, sum_x3];...
      [sum_x2, sum_x3, sum_x4] ];
RHS = [sum_S, sum_xS, sum_x2S].';
b = A\RHS;
b0 = b(1);
b1 = b(2);
b2 = b(3);

% Fit parabola to data
y = b0 + b1*x + b2*x.^2;

% Calculate R^2 coefficient of determination
R2_quad = sum( (y - mean(S)).^2 )./sum( (S - mean(S)).^2 );
disp(['R^2 = ' num2str(R2_quad)]);

% Plots results
figure(1);
set(gcf,'color','w');
scatter(x,S); hold on;
plot(x,y,'r','LineWidth',3); hold off;
legend('Raw Data','Fit-line'); legend boxoff;
legend('Location','northwest');
title(['b_0 = ' num2str(round(b0,2)) ',  b_1 = ' num2str(round(b1,2)) ...
    ',  b_2 = ' num2str(round(b2,2))]);
xlabel('x'); 
ylabel('y');
axis tight;
set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%