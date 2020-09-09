%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An Euler solver for a system of N-coupled ordinary differential equations
% (initial-value problems). See the code below for an example on how to
% implement the function.
%
% System: x'(t) = f(x(t),t) where x is a vector of N-solutions
% corresponding to N-ODEs.
%
% Inputs:
% -> func_array = a function_handle array of the form @(t,x1,x2,...,xN)
% where each array element is a function handle of @(t,x1,x2,...,xN)
% -> initial_conds = an array of numerical values for the initial
% conditions of (x1,x2,...,xN)
% -> T = total time of simulation (t is in [0,T])
%
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sols = Euler_system(func_array,initial_conds,T)
% Example: (run code without inputs)
if nargin == 0
    % 3 coupled ODEs
    f1 = @(t,x1,x2,x3) x1.^2 + cos(t);
    f2 = @(t,x1,x2,x3) x2 - x1;
    f3 = @(t,x1,x2,x3) x3 - exp(x2) - sin(x2.*x1);
    func_array = @(t,x1,x2,x3) [f1(t,x1,x2,x3); f2(t,x1,x2,x3); ...
        f3(t,x1,x2,x3)];
    initial_conds = [1.2, 1.01, 0.56];
    T = 0.5;
    plot_1 = 1;
end

x1(1) = initial_conds(1);
x2(1) = initial_conds(2);
x3(1) = initial_conds(3);

Nt = 500;
tv = linspace(0,T,Nt);
dt = tv(2)-tv(1);

% Main solver
for n = 1:length(tv)-1
    previous_step = [x1(n); x2(n); x3(n)];
    t0 = tv(n);
    initial_conds = [t0, x1(n), x2(n), x3(n)];
    forward_step = EulerStep(func_array,initial_conds,previous_step,dt);
    x1(n+1) = forward_step(1);
    x2(n+1) = forward_step(2);
    x3(n+1) = forward_step(3);
end
sols = [x1; x2; x3];

if exist('plot_1')
    close all;
    figure(1);
    set(gcf,'color','w');
    plot(tv,[x1; x2; x3],'LineWidth',3);
    xlabel('t'); ylabel('x(t)');
    legend('x_1(t)','x_2(t)','x_3(t)'); legend boxoff; 
    legend('Location','northwest');
    set(gca,'FontSize',20);
end

% Auxiliary function for forward step in time
    function forward = EulerStep(func_array,initial_conds,previous_step,dt)
        conds = num2cell(initial_conds);
        forward = previous_step + dt.*func_array(conds{:});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%