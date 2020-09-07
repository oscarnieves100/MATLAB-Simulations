%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion in 1 dimension (space-time) using the Finite-Difference Method,
% otherwise known as the solution to the Fokker-Planck Equation with zero
% drift.
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Inputs
Nx = 100; % x-grid size
xmin = -2; % left-boundary
xmax = 2; % right-boundary
x0 = 0; % initial position of particle
T = 0.25; % maximum time
D = 1.5; % diffusion coefficient
A = 4; % height of function at t=0
dx = (xmax-xmin)/Nx;
dt = 0.9*(dx^2)/(2*D); % stability condition
r = D*dt/(dx^2); % recursion constant
if r >= 1/2
    disp(['Warning: step-size dt is unstable, r = ' num2str(r)]);
    return;
end
Nt = ceil(T/dt); % time-grid size
xv = xmin:dx:xmax; % x-vector
tv = linspace(0,T,Nt); % time-vector
Lt = length(tv); % length of time
Lx = length(xv); % length of space
[TV,XV] = meshgrid(tv,xv); %grid values

% Initial conditions
p = zeros(Lx,Lt);
sigma = 0.068;
p0 = @(x) A*1/(sqrt(2*pi*sigma^2))*exp(-(x-x0).^2/(2*sigma^2));
p(:,1) = p0(xv)./(dx*trapz(p0(xv)));

%% Main loop
% Derivative matrices
M = sparse(diag(ones(Lx,1)) + r*(-2*diag(ones(Lx,1)) + ...
    diag(ones(Lx-1,1),1) + diag(ones(Lx-1,1),-1)));

for n = 1:Lt-1
    p(:,n+1) = M*p(:,n);
    
    % Normalize
    Normal = dx*trapz(p(:,n+1));
    p(:,n+1) = p(:,n+1)/Normal;
end

%% Generate plots
figure(1);
set(gcf,'color','w');
for n = 1:Lt
    plot(xv,p(:,n),'LineWidth',3);
    xlabel('x'); ylabel('p(x,t)');
    ylim([0 max(p(:,1))]);
    set(gca,'FontSize',16);
    drawnow;
end

figure(2);
set(gcf,'color','w');
surf(XV,TV,p); shading interp; grid off;
xlabel('x'); ylabel('t');
xlim([min(xv) max(xv)]); ylim([0 max(tv)]);
zlabel('p(x,t)');
set(gca,'FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%