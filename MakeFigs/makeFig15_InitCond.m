%% MAKE FIG %%

clear all;
clc;


set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

addpath ../SBPSAT
addpath ../SeismicAirgunCode

addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/598ci/FarField'
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/50ci/FarField'


%% Run Code %%

nx = 50; % number of points that air gun is discretized into per 1 m

airgunPressure = 1020; % [psi]
airgunLength = 1.2; % [m]
airgunPortArea = 12.5; % [in2]
airgunDepth = 10; % [m]
airgunVolume = airgunLength*39.3701 * airgunPortArea; % [in3]

r = 75;
rho_infty = 1e3;
c_infty = 1482; 
 
%initBubbleVol = [200 300 600 900 1200 1500 2000 2500 3000 3500 4000];
initBubbleVol_plot = [600 1200 2500 4000];
initBubbleVol = initBubbleVol_plot;

peakPressure = zeros(size(initBubbleVol));
slope = zeros(size(initBubbleVol));
riseTime = zeros(size(initBubbleVol));

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 700]);

%% DATA %%

%% Data %%

dataStr = {'219_1000cm_1030psi_598ci_DHA.csv'};
%dataStr = {'188_0750cm_1030psi_598ci_DHA.csv'};

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
 
%% Euler Air Gun Model %%

c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

tSave = [];
pSave = [];

for i = 1:length(initBubbleVol)
    disp_str = strcat('Initial Bubble Volume=',num2str(initBubbleVol(i)),'in3');
    disp(disp_str);
     
    sol = runEulerCode_initBubbleVol(nx, initBubbleVol(i),airgunPressure, airgunLength, airgunPortArea, airgunDepth);
    
    t = sol.x;
    R = sol.y(1,:);
    U = sol.y(2,:);
    m = sol.y(3,:);
    [~, solDY] = deval(sol, t);
    A = solDY(2,:);
    
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    
    [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, ... %ghost
        r + 2*airgunDepth);
    dt = 1e-5;
    tInterp = min(tDir):dt:max(tDir);
    pDirInterp = pchip(tDir, pDir, tInterp);
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp;
    
    tSave = [tSave, tInterp'];
    pSave = [pSave, pPres'];
    
    pmax(i) = max(pPres*1e-5*r);
    
    [peakPressure(i), idx] = max(pDir*1e-5*r);
    slope(i) = peakPressure(i)/((tDir(idx)-r/c_infty-tDir(1))*1000);
    riseTime(i) = (tDir(idx)-r/c_infty-tDir(1))*1000;
    
    if ismember(initBubbleVol(i), initBubbleVol_plot)
        
        
        % pressure perturbation
        subplot(3,1,[1 2]);
        plot((tInterp-r/c_infty)*1000+2, pPres*1e-5*r);
        %plot((tInterp-r/c_infty)*1000, pPres*1e-5*r);
        hold on;
        ylim([-1.8 3.8])
        
    end
    
end

%%

% rise time, slope and peak pressure
subplot(3,1,3);
plot(initBubbleVol, riseTime,'k-');
hold on;
plot(initBubbleVol, slope,'k--');
plot(initBubbleVol, peakPressure,'k:');

xlabel('Initial bubble volume (in^3)');
xlim([min(initBubbleVol) max(initBubbleVol)]);
legend('Rise Time (ms)','Slope (bar m / ms)','Peak Pressure (bar m)');

subplot(3,1,[1 2]);
ylabel('bar m');
xlabel('Time (ms)');
xlim([0 300]);

% Plot inset zoom in view of initial peak
handaxes1 = axes('Position',[0.45 0.7 0.43 0.21]);
plot(tdata*1000-tshift,pData*1e-5*r,'k');
xlim([0 12]);
ylim([0 2.5]);
xlabel('Time (ms)');
ylabel('bar m');

hold on;

k = 1;
for j = 1:length(initBubbleVol)
    if ismember(initBubbleVol(j), initBubbleVol_plot)
        plot((tSave(:,j)-r/c_inf)*1000+1,pSave(:,j)*1e-5*r,'Color',cmap(k,:)); 
        hold on;
        k = k+1;
    end
end



