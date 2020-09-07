%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion with drift (1D in space-time) using a split-step integration
% algorithm to solve the Fokker-Planck equation.
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% Inputs
Nx = 200; % x-grid size
xmin = -2; % left-boundary
xmax = 2; % right-boundary
x0 = 0; % initial position of particle
T = 3; % maximum time
mu = 3.75;
sigma1 = 0.5;
D = sigma1^2/2;
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
p(:,1) = p0(xv)./(dx.*trapz(p0(xv))); % Normalize initial condition

% Propagator matrix
diff2 = 1/(dx^2)*sparse(diag(-2*ones(Lx,1),0) + diag(ones(Lx-1,1),1) + ...
    diag(ones(Lx-1,1),-1));
diff1 = 1/(2*dx)*sparse(diag(ones(Lx-1,1),1) - diag(ones(Lx-1,1),-1));
L_op = expm( dt*(-mu*diff1 + D*diff2) );

%% Main Loop
for tt = 1:Lt-1
    % Compute forward value in time for all space
    p(:,tt+1) = L_op*p(:,tt);
    
    % Normalize probability so that total area is 1
    p(:,tt+1) = p(:,tt+1)./(dx.*trapz(p(:,tt+1)));
end

%% Generate animated plot
figure(1)
set(gcf,'color','w');
for tt = 1:Lt
    plot(xv,p(:,tt),'k','LineWidth',3);
    xlabel('x');
    ylabel('p(x,t)');
    ylim([min(p(:,1)) max(p(:,1))]);
    set(gca,'FontSize',16);
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%