% -------------------------------------------------------------------------
% This script calculates the average Fourier transform of a laser signal
% field x(t) subject to phase noise phi(t) in two different cases:
%
% (1) --> Sinusoidal field x(t) = exp(-i*w0*t)*exp(i*phi(t))
% (2) --> Gaussian field x(t) = exp(-a*t^2)*exp(-i*w0*t)*exp(i*phi(t))
%
% where w0 is the centre frequency of the laser field, and the parameter a
% which is related to the Gaussian FWHM tau via a = 4*ln(2)/tau^2.
%
% The phase noise function phi(t) is defined as a random walk (Brownian
% motion) in t, according to phi(t) = sqrt(dOm)*W(t), where dOm = 2*pi*dNu
% and dNu is the frequency linewidth of the laser source, while W(t) is a
% Wiener process defined with zero mean and variance t. In the case where
% time is taken from [-Inf,+Inf], the variance of W(t) must be defined as
% |t|. W(t) may be sampled from a Normal distribution using the definition
% W(t) = sqrt(t)*N(0,1), although in a discrete sense it is more convenient
% to define it as a sum of successive steps: 
% --> W(t+dt) = W(t) + sqrt(dt)*N(0,1)
%
% To get <X(om)>, we need <x(t)> first. Note that 
% --> < exp(i*phi(t)) > = exp(-1/2*dOm*|t|)
% which results from the calculation of the Moment-Generating function of
% W(t).
%
% It should be noted that while the Fourier Transform (FT) of (1) is a
% Lorentzian profile, the FT of (2) is a Voigt profile: a convolution
% between a Gaussian and a Lorentzian. Both profiles are subject to
% spectral broadening controlled by the laser linewidth dOm.
%
% Author: Oscar A. Nieves
% Last updated: 9/02/2021
% -------------------------------------------------------------------------
close all; clear all; clc;

%% --- INPUTS --- %%
% User-defined parameters
T = 10e-9;    % Time duraction of signal (s)
f0 = 200e12;  % Laser center frequency (Hz)
A = 1;        % Laser amplitude
FWHM = T/5;   % Full-Width at Half Maximum of Gaussian pulse
Nt = 1001;    % Sample points in time
Nf = 1001;    % Sample points in frequency
dNu = [30e5,90e5,30e6,90e6];  % Laser linewidths (Hz)

% Calculated parameters
f_end_1 = max(dNu); % Frequency around f0 for plotting
f_end_2 = 20*max(dNu);
w0 = 2*pi*f0;
w_end_1 = 2*pi*f_end_1;
w_end_2 = 2*pi*f_end_2;
t = linspace(0,T,Nt);
w1 = linspace(w0-w_end_1,w0+w_end_1,Nf);
w2 = linspace(w0-w_end_2,w0+w_end_2,Nf);
a = 4*log(2)/FWHM^2;

%% --- MAIN PROGRAM --- %%
FS = '\fontname{Palatino} ';
for n = 1:length(dNu)
    dOm = 2*pi*dNu(n);
    
    % Case 1: sinusoidal x(t)
    x1_avg = @(t) A*exp(-1i*w0*t).*exp(-dOm*abs(t)/2);
    X1_avg = @(w) sqrt(2*pi)*L(w-w0,dOm/2);
    
    % Case 2: Gaussian x(t)
    x2_avg = @(t) A*exp(-a*t.^2).*exp(-1i*w0*t).*exp(-dOm*abs(t)/2);
    X2_avg = @(w) 2*pi*V(w-w0,sqrt(2*a),dOm/2);
    
    % Normalized functions
    X1_plot(:,n) = X1_avg(w1)./max(X1_avg(w1));
    X2_plot(:,n) = X2_avg(w2)./max(X2_avg(w2));
    
    lgd{n} = [FS '\Delta\nu = ' num2str(round(dNu(n)/1e6,2)) ' MHz'];
end

%% --- PLOTS --- %%
close all;

fontS = 18;
colors = {'red','blue','green','black'};
LW = 2;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');

% Linear scale plots
subplot(221);
for n = 1:length(dNu)
    plot(w1/2/pi/1e12,X1_plot(:,n),'LineWidth',LW,'color',colors{n}); hold on;
end
hold off;
legend(lgd); legend boxoff;
axis tight; grid on;
xlabel([FS 'f (GHz)']);
ylabel([FS '|\langle' 'X(\omega)\rangle| (normalized units)']);
title([FS 'Sinusoidal signal']);
set(gca,'FontSize',fontS);

subplot(222);
for n = 1:length(dNu)
    plot(w2/2/pi/1e12,X2_plot(:,n),'LineWidth',LW,'color',colors{n}); hold on;
end
hold off;
legend(lgd); legend boxoff;
axis tight; grid on;
xlabel([FS 'f (GHz)']);
ylabel([FS '|\langle' 'X(\omega)\rangle| (normalized units)']);
title([FS 'Gaussian signal']);
set(gca,'FontSize',fontS);

% Decibel scale plots
subplot(223);
for n = 1:length(dNu)
    plot(w1/2/pi/1e12,10*log10(X1_plot(:,n)),'LineWidth',LW,'color',colors{n}); hold on;
end
hold off;
axis tight; grid on;
xlabel([FS 'f (GHz)']);
ylabel([FS '|\langle' 'X(\omega)\rangle| (normalized units)']);
title([FS 'Sinusoidal signal (dB)']);
set(gca,'FontSize',fontS);

subplot(224);
for n = 1:length(dNu)
    plot(w2/2/pi/1e12,10*log10(X2_plot(:,n)),'LineWidth',LW,'color',colors{n}); hold on;
end
hold off;
axis tight; grid on;
xlabel([FS 'f (GHz)']);
ylabel([FS '|\langle' 'X(\omega)\rangle| (normalized units)']);
title([FS 'Gaussian signal (dB)']);
set(gca,'FontSize',fontS);

%% --- FUNCTIONS --- %%
% Gaussian function
function output = G(x,sigma)
    output = 1/sqrt(2*pi*sigma^2).*exp(-x.^2/2/sigma^2);
end

% Lorentzian function
function output = L(x,gamma)
    output = 1/pi.*gamma./(gamma^2 + x.^2);
end

% Voigt profile (using Pseudo-Voigt approximation)
function output = V(x,sigma,gamma)
    fL = 2*gamma;
    fG = 2*sigma*sqrt(2*log(2));
    f = ( fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + ...
        4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5 ).^(1/5);
    eta = 1.36603*(fL/f) - 0.47719*(fL/f)^2 + 0.11116*(fL/f)^3;
    output = eta*L(x,f) + (1-eta)*G(x,f);
end