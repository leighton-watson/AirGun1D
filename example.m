%% EXAMPLE %%
%
% Example script file that shows how to run the 1D source simulations and
% plot the outputs

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

% add code directories
addpath SBPSAT/
addpath SeismicAirgunCode/

% define discretization
nx = 50; % number of grid points per 1m of source length

% source  properties
src_length = 1.2; % source length (m)
src_area = 12.5; % cross-sectional area of source (in^2)
src_depth = 7.5; % source depth (m)
src_pressure = 1000; % source pressure (psi)

% run code
tic
sol = runEulerCode(nx, src_pressure, src_length, src_area, src_depth);
toc

% bubble radius 
t = sol.x; % time vector
R = sol.y(1,:); % bubble radius

figure(1); clf;
plot(t, R);
xlabel('Time (s)');
ylabel('Bubble Radius (m)');

% acoustic pressure 
r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

U = sol.y(2,:); % bubble wall velocity (m/s)
[~,solDY] = deval(sol, t);
A = solDY(2,:); % bubble wall acceleration (m/s^2)
[tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation (Pa)
pDirBarM = pDir*1e-5*r; % convert pressure to bar m

figure(2); clf;
plot(tDir, pDirBarM);
xlabel('Time (s)');
ylabel('Pressure (bar m)');