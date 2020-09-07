%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric Brownian Motion in 1D (space-time) simulated by solving the
% corresponding Fokker-Planck equation for the probability density function
% p(x,t) using a split-step integration algorithm.
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Inputs
Nx = 500; % x-grid size
xmin = 0; % left-boundary
xmax = 20; % right-boundary
x0 = 3; % initial position of particle
T = 1.3; % maximum time
alpha = 1.5;
beta = 0.2;
D = beta^2/2;
A = 4; % height of function at t=0
dx = (xmax-xmin)/Nx;
dt = dx/10;
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

% Propagator matrix
diff2 = 1/(dx^2)*sparse(diag(-2*ones(Lx,1),0) + diag(ones(Lx-1,1),1) + ...
    diag(ones(Lx-1,1),-1));
diff1 = 1/(2*dx)*sparse(diag(ones(Lx-1,1),1) - diag(ones(Lx-1,1),-1));
L_op = expm( dt*(diag((2*D-alpha)*ones(Lx,1)) + diag(xv)*(4*D-alpha)*diff1 + D*diag(xv)^2*diff2) );

%% Main Loop
% Normalize initial p(x,t)
Normal = dx/2*(p(1,1) + 2*sum(p(2:end-1,1)) + p(end,1));
p(:,1) = p(:,1)/Normal;

for n = 1:Lt-1
    p(:,n+1) = L_op*p(:,n);
    
    % Normalize
    Normal = dx*trapz(p(:,n+1));
    p(:,n+1) = p(:,n+1)/Normal;
end

%% Plots
figure(1);
    set(gcf,'color','w');
for n = 1:Lt
    plot(xv,p(:,n),'k','LineWidth',3);
    xlabel('x');
    ylabel('p(x,t)');
    ylim([min(p(:,1)) max(p(:,1))]);
    set(gca,'FontSize',16);
    drawnow;
end

figure(2);
set(gcf,'color','w');
pcolor(TV,XV,p); shading flat; colorbar;
xlabel('t');
ylabel('x');
set(gca,'FontSize',16);

figure(3);
set(gcf,'color','w');
surf(tv,xv,p,'LineStyle','none'); grid off;
xlabel('t');
ylabel('x');
xlim([min(tv) max(tv)]);
ylim([min(xv) max(xv)]);
zlabel('p(x,t)');
set(gca,'FontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%