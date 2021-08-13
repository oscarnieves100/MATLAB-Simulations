%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solver for the SBS coupled envelope equations, as explained in
% the research paper:
%
% Oscar A. Nieves, Matthew D. Arnold, Michael J. Steel, Mikołaj K. Schmidt, 
% and Christopher G. Poulton, "Numerical simulation of noise in pulsed 
% Brillouin scattering," J. Opt. Soc. Am. B 38, 2343-2352 (2021)
%
% Inputs:
% - phase_noise/thermal_noise = either 1 or 0 (true or false) to include
%   input laser phase noise and waveguide thermal noise into the simulation
% - L = waveguide length (meters)
% - fwhm_p & fwhm_s = full-width af half-maximum of pump and Stokes pulses
%   respectively (assuming Gaussian pulses, in seconds)
% - va = acoustic group velocity (m/s)
% - n = refractive index (core waveguide material)
% - Nz = grid size in space (temporal grid size is calculated based on Nz)
% - Pp0 = input peak pump power (Watts)
% - Ps0 = input peak Stokes power (Watts)
% - tau_a = acoustic lifetime (seconds)
% - nu_laser = input laser linewidth (for pump and Stokes in Hertz)
% - lambda = input laser light central wavelength (meters)
% - f_ph = central acoustic frequency (Hertz)
% - chirp1 & chirp2 = pump and Stokes pulse constant chirp (GHz/ns)
% - alpha_dBpcm = optical loss in dB/cm
% - T = waveguide temperature (Kelvin)
% - g0 = SBS gain parameter (per meter per Watt)
% - Q2 = overlap integral (second per meter per sqrt(Watt))
% 
% Outputs:
% All outputs represent complex-valued arrays in space z and time t (rows
% correspond to z positions and columns correspond to points in time t).
% You may retrieve quantities like Stokes power by taking abs(a2).^2
% - a1 = pump envelope field (units of sqrt(Watts))
% - a2 = Stokes envelope field (units of sqrt(Watts))
% - b = acoustic envelope field (units of sqrt(Watts))
% - zv = space array (meters)
% - tv = time array (seconds)
%
% Author: Oscar Andres Nieves
% Institution: University of Technology Sydney
% Last updated: 07/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a1,a2,b,zv,tv] = SBS_solver(phase_noise,thermal_noise,L,fwhm_p,...
    fwhm_s,va,n,Nz,Pp0,Ps0,tau_a,nu_laser,lambda,...
    f_ph,chirp1,chirp2,alpha_dBpcm,T,g0,Q2)

% --- Fixed constants --- %
kB = 1.38e-23; % Boltzmann constant
h = 6.626e-34; % Planck constant
hbar = h/2/pi; % reduced Planck constant
c = 3e8;       % speed of light vacuum (m/s)

% --- Default values --- %
if ~exist('phase_noise','var') || isempty(phase_noise)
    phase_noise = 0;
end
if ~exist('thermal_noise','var') || isempty(thermal_noise)
    phase_noise = 0;
end

% --- Calculated parameters --- %
switch thermal_noise % Switch T = 0 in case of no thermal noise
    case 0
        T = 0;
    case 1
end
Gamma = 1/tau_a;       % Acoustic phonon decay rate (1/s)
v = c/n;               % optical speed (group velocity)
Om = 2*pi*f_ph;        % phonon angular freq (rad/s)
alpha_ac = Gamma/va;   % acoustic loss (1/m)
sigma = kB*T*alpha_ac; % Noise strength (J/m)
nu = c./lambda;        % Optical frequency (Hz)
w1 = 2*pi*nu;          % Pump angular frequency (rad/s)
w2 = w1 - Om;          % Stokes angular frequency (rad/s)
chirp1 = chirp1.*(fwhm_p.*1e9).^2/0.44;     % Chirp parameter data
chirp2 = chirp2.*(fwhm_s.*1e9).^2/0.44;     % Chirp parameter read/write
alpha = alpha_dBpcm*100/(10*log10(exp(1))); % Optical loss (m^{-1})
a = alpha;
tau_p = fwhm_p/2/sqrt(log(2)); % 1/e intensity at half width
tau_s = fwhm_s/2/sqrt(log(2)); % 1/e intensity at half width

% --- Gain and Overlap integrals --- %
if ~exist('Q2','var') || isempty(Q2)
    Q2 = sqrt(g0*Gamma/4/va/w2/Om);  % Stokes
end
Q1 = -conj(Q2);                       % Pump
Qa = conj(Q2);                             % Acoustic
g2 = g0;                             % Gain Stokes (energy gain)
g1 = w1/w2*conj(g2);                 % Gain Pump (energy loss)

% --- Simulation parameters --- %
transit_time = L/v;              % transit time in waveguide
d1 = 2*max([fwhm_s, fwhm_p]);    % delay of 1st Stokes pulse
simulation_time = L/v + 2*d1;
d2 = d1;

% Input Pulses (pump and Stokes)
Ap_in = @(t,C,d) sqrt(Pp0).*exp(-(1+1i*C)/2*((t - d)/tau_p).^2);
As_in = @(t,C,d) sqrt(Ps0).*exp(-(1+1i*C)/2*((t - d)/tau_s).^2);

% --- Grid Parameters --- %
zv = linspace(0,L,Nz);
dz = zv(2)-zv(1);
dt = dz/v;
Nt = ceil(simulation_time/dt);
tv = linspace(0,simulation_time,Nt);
dt = tv(2)-tv(1);
[TV,ZV] = meshgrid(tv,zv);

% --- Main Loop --- %
a1 = Ap_in(TV-ZV/v,chirp1,d1);
a2 = As_in(TV-(L-ZV)/v,chirp2,d2);
b = zeros(Nz,Nt);
ph1 = zeros(1,Nt);
ph2 = zeros(1,Nt);

% Initialize thermal field D(z,t) at t << 0
var0 = va*sigma/alpha_ac/dz;
D = sqrt(var0/2).*(randn(Nz,Nt) + 1i*randn(Nz,Nt));
b(:,1) = sqrt(dz).*D(:,1);

% Boundary Conditions (with phase noise)
switch phase_noise
    case 1
        for n = 1:Nt-1 % phase noise in laser inputs
            ph1(n+1) = ph1(n) + sqrt(2*pi*nu_laser)*sqrt(dt).*randn;
            ph2(n+1) = ph2(n) + sqrt(2*pi*nu_laser)*sqrt(dt).*randn;
        end
    case 0
        ph1 = zeros(size(tv));
        ph2 = ph1;
end
a1(1,:) = Ap_in(tv,chirp1,d1).*exp(1i*ph1);
a2(Nz,:) = As_in(tv,chirp2,d2).*exp(1i*ph2);

% Initialize overlap integrals
IG = @(a1s,a2s,nv) conj(a1s).*a2s.*exp(Gamma/2*nv*dt);
I_12 = zeros(Nz,Nt);
I_12_sum = 0;

% Time iteration
for n = 1:Nt-1
    % --- Thermal noise D --- %
    D(:,n+1) = RStep(D(:,n),dz,dt,va,Gamma,sigma);
    
    % Drift optical fields
    a1(2:Nz,n+1) = a1(1:Nz-1,n); % shift A1 to the right in z
    a2(1:Nz-1,n+1) = a2(2:Nz,n); % shift A2 to the left in z
    
    % Update Boundary Conditions
    a1(1,n+1) = Ap_in(tv(n+1),chirp1,d1).*exp(1i*ph1(n+1));
    a2(end,n+1) = As_in(tv(n+1),chirp2,d2).*exp(1i*ph2(n+1));
    
    % Temporarily store values for time evolution equations
    a1_temp = a1(:,n+1);
    a2_temp = a2(:,n+1);
    
    % Interaction integral I_12(z,t)
    if n == 1
        a1s = a1(:,1);
        a2s = a2(:,1);
    else
        a1s = [a1(:,n-1) a1_temp]; % update shifted values
        a2s = [a2(:,n-1) a2_temp]; % update shifted values
        I_12_sum = I_12_sum + IG(a1s(:,1),a2s(:,1),n-1) + ...
            IG(a1s(:,2),a2s(:,2),n);
        I_12(:,n) = dt/2*exp(-Gamma/2*n*dt)*I_12_sum;
    end
    
    a1(:,n+1) = (1 - v*a*dt/2).*a1_temp - ...
        v*dt*(g1*Gamma*conj(I_12(:,n))/4 - 1i*w1*Q1*conj(D(:,n))).*a2_temp;
    a2(:,n+1) = (1 - v*a*dt/2).*a2_temp + ...
        v*dt*(g2*Gamma*I_12(:,n)/4 - 1i*w2*Q2*D(:,n)).*a1_temp;
    b(:,n) = 1i*va*Om*Qa*I_12(:,n) + sqrt(dz).*D(:,n);
end
b(:,Nt) = 1i*va*Om*Qa*I_12(:,Nt) + sqrt(dz).*D(:,Nt);

    % --- Local functions --- %
    % Random Walk Step (Thermal noise function D(z,t))
    function [walk1] = RStep(walk0,dz,dt,va,Gamma,sigma)
        Nz0 = length(walk0);
        R1 = randn(Nz0,1);
        R2 = randn(Nz0,1);
        theta = 0.5*Gamma;
        sigma1 = va*sqrt(sigma/dz);
        R0 = 1/sqrt(2).*(R1 + 1i*R2);
        walk1 = exp(-theta*dt).*walk0 + ...
            sigma1.*sqrt( (1-exp(-2*theta*dt))/2/theta ).*R0;
    end

    % Heaviside step-function (numeric)
    function output = step_func(x)
        output = ones(size(x));
        output(x <= 0) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%