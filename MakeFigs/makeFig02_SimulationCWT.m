%% MAKE FIG 7 Pressure Simulation %%
%
% Make figure for Euler air gun "Geophysics" paper
%
% Plot mass flow rate and pressure perturbation for a range of air gun
% pressures. Plot slope, rise time and peak pressure vs initial air gun
% pressure

clear all;
clc;
%close all;

addpath ../SBPSAT
addpath ../SeismicAirgunCode

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

%figHand1 = figure(1); clf;
%set(figHand1,'Position',[100 100 600 600]);

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
%pPres = pDirInterp;

% figure(10); clf;
% plot(tInterp, pPres);
% 
% figure(11); clf;
% spectrogram(pPres,kaiser(128,18),120,128,200,'yaxis');
% xlim([0 100]/1000)
% return
    
%% plotting
figure(1); clf;
subplot(3,1,1);
plot(t*1000, R);
xmax = 500;
xlim([0 xmax])

subplot(3,1,2);
plot(t*1000, V);
xlim([0 xmax])

subplot(3,1,3);
%plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
xlim([0 xmax])

% CWT %%
Fs = 1/dt;
[wt,f] = cwt(pPres,Fs);
[m,n] = size(wt);
figure(2); clf;
h = pcolor(repmat(tInterp,m,1)*1000, repmat(f,1,n), abs(wt));
shading interp
colormap parula
ylim([5 220]);
set(h.Parent,'XTick',[0 100 200 300 400 500]);
xlim([0 xmax]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%% other plots %%

V = 4/3*pi*R.^3;
Vdot = 4*pi*R.^2.*U;
Vddot = 8*pi*R.*U.^2 + 4*pi*R.^2.*A;


figHand3 = figure(3); clf;
set(figHand3,'Position',[100 100 600 900]);
subplot(4,1,1);
plot(t*1000,V);
xlim([0 xmax]);
hold on;
xlabel('Time, t (ms)');
ylabel('Volume (m^3)');
ylim([0 1.5])
% h = text(0.006*1000, 1.25, '(a)');
% set(h,'FontSize',24);
% set(h,'FontWeight','bold');
%plot(t*1000,mass);

% subplot(4,1,2);
% plotyy(t*1000,Vdot,t*1000,Vddot);
% xlim([0 xmax]);

subplot(4,1,2);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
%plot((tInterp)*1000, pPres*1e-5*r);
xlim([0 xmax]);
set(h.Parent,'YTick',[-4 0 4]);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
ylim([-4 6])
% h = text(0.006*1000, 4.6, '(b)');
% set(h,'FontSize',24);
% set(h,'FontWeight','bold');

subplot(4,1,[3 4]);
h = pcolor(repmat(tInterp-r/c_inf,m,1)*1000, repmat(f,1,n), abs(wt));
%h = pcolor(repmat(tInterp,m,1)*1000, repmat(f,1,n), abs(wt));
shading interp
colormap parula
ylim([5 220]);
%set(h.Parent,'XTick',[0 50 100 150 200]);
xlim([0 xmax]);
xlabel('Time, t-r/c_\infty (ms)');
ylabel('Frequency (Hz)');
% h = text(0.006*1000, 205, '(c)');
% set(h,'FontSize',24);
% set(h,'FontWeight','bold');
% set(h,'Color',[1 1 1]);
%keyboard;


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
% h = text(1.7, 5.4, '(d)');
% set(h,'FontSize',24);
% set(h,'FontWeight','bold');
%set(h,'Color',[1 1 1]);

h = text(0.06, 5, 'peak pressure');
set(h,'FontSize',24);

h = text(1.1, 1.5, 'rise');
set(h,'FontSize',24);
h = text(1.1, 0.7, 'time');
set(h,'FontSize',24);

