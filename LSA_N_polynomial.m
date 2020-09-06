%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares-Analysis (LSA) for n-order polynomial fit
%
% Description: generates a random data-set around a polynomial curve with
% preset values and then uses LSA matrix methods to find the coefficients
% of the line of best-fit, and overlays the fit-line to the dataset:
% Equation --> y = b0 + b1*x + b2*x^2 + b3*x^3 + ... + b{n}*x^{n}
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Generate dataset
x = 0:0.01:2; % x-axis

% Set 'guess' for order of polynomial needed
N_guess = 7;

% Generate N random coefficients for polynomial using a Uniform Distribution
% between the values [a,b]
a = -1;
b = 1;
N = 7; % order of dataset polynomial
coeffs = a + (b - a)*rand(N+1,1);
coeffs_ran = a + (b - a)*rand(N+1,1);

% dataset
S = 0;
for k = 1:N+1
   S = S + coeffs(k)*x.^(k-1) + coeffs_ran(k).*randn(size(x));  
end

%% Polynomial LSA
n = N_guess;
for k = 1:2*n+1
   row_A(k) = sum( x.^(k-1) );
end

for k = 1:n+1
    A(k,:) = row_A(k:end-n+k-1);
    RHS(k) = sum( S.*x.^(k-1) );
end

% convert b into a column vector
if size(RHS,1) == 1
    RHS = RHS.';
end

% Compute coefficient vector
b = A\RHS;

% Fit polynomial to data
y = 0;
for k = 1:n+1
   y = y + b(k).*x.^(k-1);
end

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
xlabel('x'); 
ylabel('y');
axis tight;
set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%