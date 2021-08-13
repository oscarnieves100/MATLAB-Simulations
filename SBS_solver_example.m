%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solver for the SBS coupled envelope equations, as explained in
% the research paper:
%
% Oscar A. Nieves, Matthew D. Arnold, Michael J. Steel, Mikołaj K. Schmidt, 
% and Christopher G. Poulton, "Numerical simulation of noise in pulsed 
% Brillouin scattering," J. Opt. Soc. Am. B 38, 2343-2352 (2021)
%
% Example of how to use the SBS_solver() function.
%
% * Important note: you must set the random number generator seed (here
% denoted by the variabled seed0) before doing multiple Monte Carlo runs.
% This is to ensure consistency between computations.
%
% Author: Oscar Andres Nieves
% Institution: University of Technology Sydney
% Last updated: 07/07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% --- Fixed Consants --- %%
% Noise properties (1 = Include noise, 0 = Don't include noise)
phase_noise = 1;    % Input laser phase noise
thermal_noise = 1;  % Thermal noise inside waveguide

%% --- Input parameters --- %%
% Random number generator seed
seed0 = 1;
rng(seed0);

% Waveguide properties
T = 300;              % Temperature: set = 0 if no thermal noise is wanted
alpha_dBpcm = 0.05;   % Optical loss (dB/cm)
n = 2.44;             % Refractive index
L = 0.5;              % Length (m)
tau_a = 5.3e-9;       % acoustic lifetime (s)
va = 2500;            % acoustic speed (group velocity m/s)
g0 = 423;             % SBS gain (1/m/W)

% Input laser properties
Pp0 = 1;              % Pump peak power (W) (Write pulse)
Ps0 = 1e-6;           % Stokes peak power (W)
lambda = 1550e-9;     % optical wavelength (m)
f_ph = 7.8e9;         % phonon freq (Hz)
nu_laser = 100e3;     % Laser linewidth (Hz)
chirp1 = 0;           % Chirp rate input pump (GHz/ns)
chirp2 = 0;           % Chirp rate input Stokes (GHz/ns)
fwhm_p = 2e-9;        % Pump FWHM (s)
fwhm_s = 1e-9;        % Stokes FWHM (s)

% Numerical simulation parameters
N = 1;                % Number of ensembles (Monte Carlo runs)
Nz = 501;             % Spatial grid size

% Animation options
anim_active = 0;      % 1 = produce animation, 0 = skip animation
title1 = 'SBS_solver_animation.mp4';   % Animation title
data_title = 'SBS_solution_data.mat';  % Data ttle

%% --- Run SBS Solver --- %%
P1avg = 0; P2avg = 0; Paavg = 0;
P1sq = 0; P2sq = 0; Pasq = 0;

% Monte Carlo runs
tic;
for kk = 1:N
    % Solve envelope fields
    [a1,a2,b,zv,tv] = SBS_solver(phase_noise,thermal_noise,L,fwhm_p,...
        fwhm_s,va,n,Nz,Pp0,Ps0,tau_a,nu_laser,lambda,...
        f_ph,chirp1,chirp2,alpha_dBpcm,T,g0);
    
    % Store single realization
    if kk == 1
        P1 = abs(a1).^2;
        P2 = abs(a2).^2;
        Pa = abs(b).^2;
    end
    
    % Calculate average values
    P1avg = P1avg + abs(a1).^2/N;
    P2avg = P2avg + abs(a2).^2/N;
    Paavg = Paavg + abs(b).^2/N;
    P1sq = P1sq + abs(a1).^4/N;
    P2sq = P2sq + abs(a2).^4/N;
    Pasq = Pasq + abs(b).^4/N;
end
comp_time = toc;
disp(['Computation time = ' num2str(comp_time/60) ' min']);

% Calculate variance and standard deviation of the powers
P1_std = sqrt( abs(P1sq - P1avg.^2) );
P2_std = sqrt( abs(P2sq - P2avg.^2) );
Pa_std = sqrt( abs(Pasq - Paavg.^2) );

%% --- Plot Results --- %%
close all;
colors = {'black','#D95319','#0072BD','#77AC30'};
Nt = length(tv);

% Plot characteristics
step = 2.*ceil(Nt/50);
fontS = 16;
FS = '\fontname{Palatino} ';
x_ticks = linspace(0,L,10);
[ZV1,TV1] = meshgrid(1e2*zv,tv(:,1:step:Nt));
[ZV2,TV2] = meshgrid(1e2.*zv,tv);
az_angle1 = 30;
el_angle1 = 50;
fontSS = 20;

% --- Figures --- %
yLIMS = [0,max(tv.*1e9)];
xLIMS = [0,L*1e2];
xticks_v = linspace(0,L*1e2,6);

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');
set(0,'DefaultAxesTitleFontWeight','normal');

AX1 = subplot(131);
pp1 = waterfall(ZV1,TV1.*1e9,P1(:,1:step:Nt).'); grid off;
pp1.EdgeColor = colors{3};
pp1.LineWidth = 2;
view(az_angle1,el_angle1);
xlim(xLIMS);
ylim(yLIMS);
xticks(xticks_v);
xlabel([FS 'z (cm)']);
ylabel([FS 't (ns)']);
zlabel([FS 'Pump power (W)']);
%title([FS 'P_1(z,t) (max ' num2str(maxP1) ' W)']);
set(gca,'FontSize',fontSS);

AX2 = subplot(132);
pp2 = waterfall(ZV1,TV1.*1e9,P2(:,1:step:Nt).'); grid off;
pp2.EdgeColor = colors{2};
pp2.LineWidth = 2;
view(az_angle1,el_angle1);
xlim(xLIMS);
ylim(yLIMS);
xticks(xticks_v);
xlabel([FS 'z (cm)']);
ylabel([FS 't (ns)']);
zlabel([FS 'Stokes power (W)']);
%title([FS 'P_2(z,t) (max ' num2str(maxP2) P2title]);
set(gca,'FontSize',fontSS);

AX3 = subplot(133);
pp3 = waterfall(ZV1,TV1.*1e9,Pa(:,1:step:Nt).'); grid off;
pp3.EdgeColor = colors{1};
view(az_angle1,el_angle1);
xlim(xLIMS);
ylim(yLIMS);
xticks(xticks_v);
xlabel([FS 'z (cm)']);
ylabel([FS 't (ns)']);
zlabel([FS 'Acoustic power (W)']);
%title([FS 'P_a(z,t) (max ' num2str(maxPa) Patitle]);
set(gca,'FontSize',fontSS);

save(data_title);

%% --- Animation --- %%
if anim_active == 1
    % Plot characteristics
    close all;
    step = 2.*ceil(Nt/50);
    fontS = 22;
    FS = '\fontname{Palatino} ';
    x_ticks = linspace(0,L,10);
    [ZV1,TV1] = meshgrid(zv,tv(:,1:step:Nt));
    az_angle = 60;
    el_angle = 60;
    fontSS = fontS;
    maxP1 = max(max(P1));
    maxP2 = max(max(P2));
    maxPa = max(max(Pa));
    medZ = ceil(median(1:length(zv)));
    disp(['max P1 = ' num2str(max(max(P1))) ' W']);
    disp(['max P2 = ' num2str(max(max(P2))) ' W']);
    disp(['max Pa = ' num2str(max(max(Pa))) ' W']);
    disp(['max P1 dev = ' num2str(max(max(P1_std))) ' W']);
    disp(['max P2 dev = ' num2str(max(max(P2_std))) ' W']);
    disp(['max Pa dev = ' num2str(max(max(Pa_std))) ' W']);
    
    % Animation
    step1 = 1;
    vid1 = VideoWriter(title1,'MPEG-4');
    vid1.Quality = 100;
    open(vid1);
    figure(1);
    set(gcf,'color','w');
    set(gcf, 'Position',  [100, 100, 1200, 400]);
    for n = 1:step1:Nt
        subplot(131);
        h1 = area(zv, P1(:,n));
        h1.LineWidth = 0.5;
        h1.FaceColor = 'blue';
        h1.FaceAlpha = 0.5;
        h1.EdgeColor = 'black';
        xlim([0 L]);
        ylim([0 maxP1]);
        xlabel([FS 'z (m)']);
        ylabel([FS 'Pump power (W)']);
        title([FS 't = ' num2str(round(1e9.*tv(n),2)) ' (ns)']);
        set(gca,'FontSize',fontS);
        
        subplot(132);
        h2 = area(zv, P2(:,n));
        h2.LineWidth = 0.5;
        h2.FaceColor = 'red';
        h2.FaceAlpha = 0.5;
        h2.EdgeColor = 'black';
        xlim([0 L]);
        ylim([0 maxP2]);
        xlabel([FS 'z (m)']);
        ylabel([FS 'Stokes power (W)']);
        set(gca,'FontSize',fontS);
        
        subplot(133);
        h3 = area(zv, Pa(:,n));
        h3.LineWidth = 0.5;
        h3.FaceColor = 'green';
        h3.FaceAlpha = 0.5;
        h3.EdgeColor = 'black';
        xlim([0 L]);
        ylim([0 maxPa]);
        xlabel([FS 'z (m)']);
        ylabel([FS 'Acoustic power (W)']);
        set(gca,'FontSize',fontS);
        
        % Export frame to video
        drawnow;
        frame = getframe(gcf);
        writeVideo(vid1,frame);
    end
    close(vid1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%