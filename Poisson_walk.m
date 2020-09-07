%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A simple Poisson random walk simulator in 1-dimension
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Parameters
T = 0.5;     % simulation time
sims = 5;    % number of walks
dt = 1e-3;   % step-size
tv = 0:dt:T; % time array
X0 = 0;      % initial position
lam = 20;   % lambda parameter
Lt = length(tv);
Y = zeros(Lt,sims); % store random walks

%% Random walk loop
for N = 1:sims
    X = zeros(Lt,1);
    X(1) = X0;
    for n = 1:Lt-1
        Dr = poissrnd(lam*dt);
        X(n+1) = X(n) +  Dr;
    end
    Y(:,N) = X;
end

%% Plots
figure(1);
set(gcf,'color','w');
plot(tv,Y,'LineWidth',3);
xlabel('Time (s)');
ylabel('X(t)');
title('Poisson walks');
set(gca,'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%