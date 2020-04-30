%% MAKE FIG 13 TURBULENT CORRECTIONS %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Compare data and simulation results for model with turbuluent corrections
% added.
%
% For information about the data see Ronen and Chelminski (2018) A next 
% generation seismic source with low frequency signal and low 
% environmental impact, 80th EAGE Conference & Exhibition. 
% doi:10.3997/2214-4609.201800745

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

% add directories
addpath ../Data
addpath ../SeismicAirgunCode/
addpath ../SBPSAT/
addpath ../sbplib/
addpath ../SBPSAT/Turbulence/


figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 700]);

%% Data %%

dataStr = {'219_1000cm_1030psi_598ci_DHA.csv'};

tshift = 90.27;
data = csvread(dataStr{1});
pData = data(:,1);
    
dt = 31.25e-6; % sampling internal [ms]
Fs = 1/dt; % sampling frequency [Hz]
k = length(data);
tdata = 0:dt:(k-1)*dt;

r = 75; % distance from source to receiver [m]

figure(1);
subplot(3,1,[1 2]);
plot(tdata*1000-tshift,pData*1e-5*r,'k');
hold on;
xlim([0 300])
ylabel('bar m');

dmax = max(pData*1e-5*r);

%% Run 1D Air Gun Simulation %%

nx = 50; % number of grid points per 1 m of air gun length

c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

in2_m2 = 0.00064516; % conversion from in^2 to m^2
psi_pa = 6894.76; % conversion from psi to pa

gamma = 1.4; % ratio of heat capacities
Q = 287.06; % specific gas constant for dry air [J/kgK]
T_inf = 288; % temperature assumed constant throughout the system [K]

aP = 1030; % air gun pressure
aL = 1.2; % air gun length [m]
aA = 12.5; % air gun port area [in^2] % cross-sectional area = port area
aD = 10; % air gun depth [m]

beta = [0 2 1 5]; % damping parameter
dissp_boolean = logical([1 1 0 0]); % constant or time varying dissipation

tSave = [];
pSave = [];

for i = 1:length(beta) % beta and C need to be the same length

    sol = runEulerCode_turbCorrec(nx, aP, aL, aA, aD, 0, dissp_boolean(i), beta(i));
        
    t = sol.x; % time
    
    % bubble properties
    R = sol.y(1,:); % bubble radius [m]
    U = sol.y(2,:); % bubble wall velocity [m/s] 
    [~,solDY] = deval(sol, t);
    A = solDY(2,:); % acceleration
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, ... %ghost
        r + 2*aD);
    dt = 1e-5; % time step for interpolation
    tInterp = min(tDir):dt:max(tDir); % interpolate onto a common time vector
    pDirInterp = pchip(tDir, pDir, tInterp); 
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp; % direct and ghost
    
    figure(1);
    subplot(3,1,[1 2]);
    plot((tInterp-r/c_inf)*1000+2,pPres*1e-5*r,'Color',cmap(i,:)); 
    hold on;
    
    pmax(i) = max(pPres*1e-5*r);
        
    subplot(3,1,3);
    plot(t*1000, U); hold on;
    
    tSave = [tSave, tInterp'];
    pSave = [pSave, pPres'];
    
end

figure(1)
subplot(3,1,3);
xlim([0 300]);
ylabel('m/s');
ylim([-15 35]);
xlabel('Time (ms)');

subplot(3,1,[1 2]);
xlim([0 300]);
ylim([-1.5 4])
ylabel('bar m');

% Plot inset zoom in view of initial peak
handaxes1 = axes('Position',[0.45 0.7 0.43 0.21]);
plot(tdata*1000-tshift,pData*1e-5*r,'k');
xlim([0 12]);
ylim([0 2.5]);
xlabel('Time (ms)');
ylabel('bar m');

hold on;

for j = 1:length(beta)
    plot((tSave(:,j)-r/c_inf)*1000+2,pSave(:,j)*1e-5*r,'Color',cmap(j,:)); 
    hold on;
end



