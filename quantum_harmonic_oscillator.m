%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantum harmonic oscillator
%
% This program solves the 1-dimensiona quantum harmonic oscillator problem
% (quadratic potential) using both the analytic solutions in terms of
% Hermite polynomials and an eigenmode solver implementing finite
% differences.
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% Inputs
m = 1;
hbar = 1;
omega = 0.5;
dx = 0.1;
x = -5:dx:5;
N = length(x);

%% Numerical solution (Eigenmode solver)
% Derivative matrix
Dup = diag(ones(N-1,1),1);
Dmid = -2*diag(ones(N,1),0);
Dlow = diag(ones(N-1,1),-1);
Dxx = 1/(dx^2).*( Dlow + Dmid + Dup );

% Linear operator A
A = -hbar^2/2/m * Dxx + 1/2*m*omega^2*diag(x.^2,0);

% Eigenvalue solver
[eigenfunc, eigenval] = eig(A);
n = [0, 1, 2, 3]; % first three energy levels
for kk = 1:length(n)
    wavefunction(:,kk) = eigenfunc(:,kk);
    Energy(kk) = eigenval(kk,kk);
end

% print energies
disp(Energy);

% Normalize wavefunctions using trapezoidal rule
for kk = 1:length(n)
    Prob = wavefunction(:,kk).*conj(wavefunction(:,kk));
    Integral_Prob = dx/2.*sum( Prob(1:end-1) + Prob(2:end) );
    wavefunction(:,kk) = wavefunction(:,kk)./sqrt( Integral_Prob );
end

%% Analytic solutions (Using Hermite Polynomials)
psi_n = @(xv,nv) (m*omega/pi/hbar).^(1/4).*1/sqrt(2.^(nv).*factorial(nv)).*...
    hermiteH(nv,sqrt(m*omega/hbar).*xv).*exp(-m*omega/2/hbar.*xv.^2);
analytic_energies = hbar*omega*(n + 1/2);
disp('analytic energies:');
disp(analytic_energies);

%% Plot eigenfunctions
figure(1);
subplot(121);
set(gcf,'color','w');
for kk = 1:length(n)
    plot(x,wavefunction(:,kk),'LineWidth',1.5); hold on;
    legend_entries{kk} = ['n = ' num2str( n(kk) )];
end
hold off;
legend(legend_entries); legend boxoff;
xlabel('x');
ylabel('\psi_n(x)');
title('Eigenmode Solver');
set(gca,'FontSize',20);

subplot(122);
markers = {'s','--','o','*'};
for kk = 1:length(n)
    plot(x,psi_n(x,n(kk)),'LineWidth',1.5); hold on;
    legend_entries{kk} = ['n = ' num2str( n(kk) )];
end
hold off;
legend(legend_entries); legend boxoff;
xlabel('x');
ylabel('\psi_n(x)');
title('Analytic Solutions');
set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%