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
% addpath sbplib/ % add path to sbplib

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

figure(1); clf; subplot(3,1,1);
plot(t, R);
xlabel('Time (s)');
ylabel('m'); xlim([0 2]);
title('(a) Bubble Radius');

% acoustic pressure 
r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

U = sol.y(2,:); % bubble wall velocity (m/s)
[~,solDY] = deval(sol, t);
A = solDY(2,:); % bubble wall acceleration (m/s^2)
[tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation (Pa)
pDirBarM = pDir*1e-5*r; % convert pressure to bar m

figure(1); subplot(3,1,2);
plot(tDir, pDirBarM);
xlabel('Time (s)');
ylabel('bar m'); xlim([0 2]);
title('(b) Acoustic Pressure');

% pressure inside source
t = sol.x; % time
x = [0:ceil(src_length*nx)]./nx; % space vector
[T,X] = meshgrid(t,x); % create mesh for space-time plots

% initialize matrices
rho = zeros(length(x), length(t));
rhov = zeros(length(x), length(t));
e = zeros(length(x), length(t));

for i = 1:length(x) % extract air gun properties
    rho(i,:) = sol.y(3*i+2,:); % density
    rhov(i,:) = sol.y(3*i+3,:); % density * velocity
    e(i,:) = sol.y(3*i+4,:); % internal energy
end

gamma = 1.4; % ratio of heat capacities
v = rhov./rho; % velocity [m/s]
p = (gamma-1)*(e-0.5*rho.*v.^2); % pressure

% pressure in source
figure(1); subplot(3,1,3);
pa2psi = 0.000145038; % conversion from pa to psi
h = surf(X,T*1000,p*pa2psi);
view(2); shading interp
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Pressure');
cb = colorbar; cb.Label.String = 'psi';
ylim([0 20]); xlim([0 src_length]);
set(h.Parent,'XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
title('(c) Source Pressure');


