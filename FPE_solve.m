%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D Fokker-Planck Equation Solver. Uses a split-step integration algorithm
% to solve the equation for a probability density p(x,t) for given input
% drift mu(x,t) and volatility sigma(x,t) functions:
%
% dp/dt = -d/dx[mu(x,t)*p] + d^2/dx^2[sigma^2(x,t)/2*p]
%
% which corresponds to the Ito SDE:
%
% dX = mu(X,t)dt + sigma(X,t)dW
% 
% This is done by casting the problem into the matrix form: 
% 
% d/dt[p(x,t)] = H(x,t)*p(x,t), which is solved via the equation
%
% p(x,t+dt) = exp(dt*H(x,t))*p(x,t)
%
% *Note: because the spatial grid is restricted to the finite domain [a,b],
% the function p(x,t) will deform near the boundaries to ensure the total
% probability p(x,t0) at any time t0 is always 1 (due to normalization
% requirements). To model realistic boundary-less systems, it is
% recommended to set a or b very large depending on the problem.
%
% Inputs:
% -> mu(X,t) and sigma(X,t) (function_handles @(x,t))
% -> p0(x) = initial density function (function_handle @(x))
% -> T = final time of simulation
% -> xlimits = limits on x-axis (e.g. [a,b])
% 
% Made by: Oscar A. Nieves
% Made in: 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = FPE_solve(mu,sigma,p0,T,xlimits)
% If no inputs, use default values (example: Geometric Brownian Motion)
if nargin == 0
    u = 1.0;
    os = 0.25;
    mu = @(x,t) u*x;
    sigma = @(x,t) os*x;
    p0 = @(x) 1/(2*pi*os.^2).*exp(-(x-1).^2/2/os^2);
    T = 1;
    xlimits = [-5,20];
end

% Create arrays
Nx = 200;
Nt = 300;
xv = linspace(xlimits(1),xlimits(2),Nx);
tv = linspace(0,T,Nt);
dx = xv(2)-xv(1);
dt = tv(2)-tv(1);
p = zeros(Nx,Nt);
p(:,1) = p0(xv)./(dx*trapz(p0(xv)));

% Finite difference matrices
Dxx = 1/(dx^2)*sparse(diag(-2*ones(Nx,1),0) + diag(ones(Nx-1,1),1) + ...
    diag(ones(Nx-1,1),-1));
Dx = 1/(2*dx)*sparse(diag(ones(Nx-1,1),1) - diag(ones(Nx-1,1),-1));

% Propagation operator
H = @(x,t) -( Dx*diag(mu(x,t)) + diag(mu(x,t))*Dx ) + ...
    1/2*( Dxx*diag( sigma(x,t).^2 ) + 2*Dx*diag( sigma(x,t).^2 )*Dx + ...
    diag( sigma(x,t).^2 )*Dxx );
L_op = @(x,t) expm( dt.*H(x,t) );

% Main Loop
for n = 1:Nt-1
    % Update value of p(x,t)
    p(:,n+1) = L_op(xv,tv(n))*p(:,n);
    
    % Normalize
    p(:,n+1) = p(:,n+1)/(dx*trapz(p(:,n+1)));
end

% Plots
[XV,TV] = meshgrid(xv,tv);

figure(1);
set(gcf,'color','w');
for n = 1:Nt
    plot(xv,p(:,n),'k','LineWidth',3);
    xlabel('x');
    ylabel('p(x,t)');
    ylim([min(p(:,1)) max(p(:,1))]);
    set(gca,'FontSize',18);
    drawnow;
end

figure(2);
set(gcf,'color','w');
pcolor(TV,XV,p.'); shading flat; colorbar;
xlabel('t'); ylabel('x');
set(gca,'FontSize',18);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%