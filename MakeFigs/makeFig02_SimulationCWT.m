%% MAKE FIG 2 SIMULATION CWT %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Display 1D air gun simulation results. Plot bubble volume time series,
% acoustic pressure time series and continuous wavelet transform

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);

% add code directories
addpath ../SBPSAT
addpath ../SeismicAirgunCode

%% Run Euler Air Gun Simulation %%

nx = 100; % number of grid points per 1 m of air gun length

tmin = 0;
tmax = 10;

r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

in2_m2 = 0.00064516; % conversion from in^2 to m^2
psi_pa = 6894.76; % conversion from psi to pa

gamma = 1.4; % ratio of heat capacities
Q = 287.06; % specific gas constant for dry air [J/kgK]
T_inf = 288; % temperature assumed constant throughout the system [K]

aP = 2000; % air gun pressure for slope/rise time plot
aL = 0.6; % air gun length [m]
aA = 16; % air gun port area [in^2] % cross-sectional area = port area
aD = 7.5; % air gun depth [m]

% run solve
sol = runEulerCode(nx, aP, aL, aA, aD);

t = sol.x; % time

% bubble properties
R = sol.y(1,:); % bubble radius [m]
V = 4/3*pi*R.^3; % bubble volume [m^3]
U = sol.y(2,:); % bubble wall velocity [m/s]
mass = sol.y(3,:); % bubble mass [kg]
[~,solDY] = deval(sol, t);
A = solDY(2,:); % acceleration
[tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
[tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r + 2*aD); %ghost
    
dt = 1e-4;
tInterp = min(tDir):dt:max(tDir);
pDirInterp = pchip(tDir, pDir, tInterp);
pGhostInterp = pchip(tGhost, pGhost, tInterp);
pPres = pDirInterp - pGhostInterp;

%%% plotting %%%    
figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 900]);

% bubble volume
V = 4/3*pi*R.^3;
subplot(4,1,1);
plot(t*1000,V);
xmax = 500;
xlim([0 xmax]);
hold on;
xlabel('Time, t (ms)');
ylabel('Volume (m^3)');
ylim([0 1.5])

% acoustic pressure
subplot(4,1,2);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
xlim([0 xmax]);
set(h.Parent,'YTick',[-4 0 4]);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
ylim([-4 6])

% CWT %%
Fs = 1/dt;
[wt,f] = cwt(pPres,Fs);
[m,n] = size(wt);
subplot(4,1,[3 4]);
h = pcolor(repmat(tInterp-r/c_inf,m,1)*1000, repmat(f,1,n), abs(wt));
shading interp
colormap parula
ylim([5 220]);
xlim([0 xmax]);
xlabel('Time, t-r/c_\infty (ms)');
ylabel('Frequency (Hz)');

handaxes1 = axes('Position',[0.47 0.26 0.4 0.2]);
set(handaxes1,'YColor',[0 0 0]);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
xlim([0 2])
ylim([0 6])
set(h.Parent,'YColor',[1 1 1]);
set(h.Parent,'XColor',[1 1 1]);
xlabel('Time (ms)');
ylabel('\Delta p (bar m)');
hold on;
[xx,idx] = max(pPres*1e-5*r);
g = hline(xx);
set(g,'Color','k');
set(g,'LineStyle','--');
g = vline((tInterp(idx)-r/c_inf)*1000);
set(g,'Color','k');
set(g,'LineStyle','--');

h = text(0.06, 5, 'peak pressure');
set(h,'FontSize',24);
h = text(1.1, 1.5, 'rise');
set(h,'FontSize',24);
h = text(1.1, 0.7, 'time');
set(h,'FontSize',24);

