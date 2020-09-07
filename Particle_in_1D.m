%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantum particle in a box (1-dimension)
%
% This program solves the particle in a box problem on a symmetric
% potential of length L using both the analytic solution and the eigenmode
% solution (e.g. using a matrix of finite differences for the Hamiltonian
% operator H in the equation: H*psi(x) = E*psi(x))
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

% Inputs (Setting Planck's constant to 1)
hbar = 1;
m = 1;
L = 4;
x = -L/2:0.01:L/2;
dx = x(2)-x(1);
N = length(x);

% Derivative matrix Dxx
main_diag = diag(-2*ones(N,1),0);
upper_diag = diag(ones(N-1,1),1);
lower_diag = diag(ones(N-1,1),-1);
Dxx = 1/dx^2 * (lower_diag + main_diag + upper_diag);

% Linear operator on the left (matrix A)
A = -hbar^2/2/m * Dxx;

% Eigenvalues and eigenvectors of A
[A_eigenvectors, A_eigenvalues] = eig(A);

% Extract energies and solutions psi(x)
sols = 3;
for n = 1:sols
    Energies(n) = A_eigenvalues(n,n);
    psi(:,n) = A_eigenvectors(:,n);
end

% Normalize eigenfunctions
for n = 1:3
    psi_sq = psi(:,n).*conj(psi(:,n));
    normalization = dx/2 * sum( psi_sq(1:end-1) + psi_sq(2:end) );
    psi(:,n) = psi(:,n)/sqrt(normalization);
end

% Plot first 3 solutions (energy levels)
figure(1);
set(gcf,'color','w');
subplot(121);
plot(x, psi(:,1:sols), 'LineWidth',3);
xlabel('x');
ylabel('\psi(x)');
legend('\psi_1(x)','\psi_2(x)','\psi_3(x)'); legend boxoff;
set(gca,'FontSize',20);

subplot(122);
plot(x, psi_analytic(x,1),x, psi_analytic(x,2),x, psi_analytic(x,3), 'LineWidth',3);
xlabel('x');
ylabel('\psi(x)');
legend('\psi_1(x)','\psi_2(x)','\psi_3(x)'); legend boxoff;
set(gca,'FontSize',20);

% Analytic solution
function psi_out = psi_analytic(xv,nv)
L = 4;
if rem(nv,2) == 0 
    psi_out = sqrt(2/L)*sin(nv*pi.*xv/L);
else
    psi_out = sqrt(2/L)*cos(nv*pi.*xv/L);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%