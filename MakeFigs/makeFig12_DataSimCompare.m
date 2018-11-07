%% DATA SIM COMPARE %%

clear all;
clc;
%close all;

addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/598ci/FarField'
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/50ci/FarField'


set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');


figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 500]);

dt = 31.25e-6; % sampling internal [ms]
Fs = 1/dt; % sampling frequency [Hz]

r = 75; % distance from source to receiver [m]

%% Vary pressure %%
% dataStr = {'155_0500cm_1030psi_598ci_DHA.csv',...
%     '156_0500cm_1040psi_598ci_DHA.csv',...
%     '157_0500cm_1030psi_598ci_DHA.csv',...
%     '158_0500cm_1020psi_598ci_DHA.csv',...
%     '159_0500cm_1030psi_598ci_DHA.csv',...
%     '160_0500cm_1030psi_598ci_DHA.csv'};

%dataStr = {'188_0750cm_1030psi_598ci_DHA.csv'};
dataStr = {'219_1000cm_1030psi_598ci_DHA.csv'};


tshift = 91-0.73;
%tshift = [1.094 0.7188 1 0.75 0.8438 0.9688]+91; %[0 1.269 1.45 0.3 0.77];    
pshift = [0]; %[-0.005 0.004 -0.004 -0.004 -0.008];
data = csvread(dataStr{1});
pData = data(:,1);

k = length(data);
tdata = 0:dt:(k-1)*dt;

plot(tdata*1000-tshift,pData*1e-5*r+pshift,'k');
hold on;
xlim([0 400])




%% Simulation %%

addpath ../SBPSAT
addpath ../SeismicAirgunCode

%% Run Euler Air Gun Simulation %%

nx = 100; % number of grid points per 1 m of air gun length

tmin = 0;
tmax = 1;

r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

in2_m2 = 0.00064516; % conversion from in^2 to m^2
psi_pa = 6894.76; % conversion from psi to pa

gamma = 1.4; % ratio of heat capacities
Q = 287.06; % specific gas constant for dry air [J/kgK]
T_inf = 288; % temperature assumed constant throughout the system [K]

aP = 1020;
aL = 1.2; % air gun length [m]
aA = 12.5; % air gun port area [in^2] % cross-sectional area = port area
aD = 10; % air gun depth [m]


sol = runEulerCode(nx, aP, aL, aA, aD);

t = sol.x; % time

% bubble properties
R = sol.y(1,:); % bubble radius [m]
U = sol.y(2,:); % bubble wall velocity [m/s]
m = sol.y(3,:); % bubble mass [kg]
[~,solDY] = deval(sol, t);
A = solDY(2,:); % acceleration
[tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
%%% include ghost % mismatch between data and simulation during decay
%%% of peak may be due to treatment of the ghost. Have not been
%%% including the ghost
[tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, ... %ghost
    r + 2*aD);
dt = 1e-5;
tInterp = min(tDir):dt:max(tDir);
pDirInterp = pchip(tDir, pDir, tInterp);
pGhostInterp = pchip(tGhost, pGhost, tInterp);
pPres = pDirInterp - pGhostInterp;

% plot((tInterp-r/c_inf)*1000+2,pPres*1e-5*r,'Color',cmap(1,:)); %,'LineStyle','--','Color',cmap(i,:));

%xlim([-0.5 8.5])
xlim([0 400]);
ylim([-2.5 4]);
xlabel('Time (ms)');
ylabel('bar m');
title('Acoustic Pressure');

%% Plot lumped paramter model %%

m_in = 39.3701; % conversion from m to in
aV = aL*m_in * aA; % air gun volume [in^3]
input = [aP, aV, aA]; % inputs for lumped parameter model [pressure, volume, port area]
physConst = physical_constants(aD,r); % save physical constants. Specific depth and distance from source to receiver
output = AirgunBubbleSolveOutput(input, physConst, false); % solve lumped parameter model

% longer time series of pressure perturbation

h = plot((output.tPres-r/physConst.c_infty)*1000+2,output.pPres*1e-5*r,...
    'Color',cmap(2,:),'LineStyle','-');
%h.Color(4) = alpha;
%xlim([tmin tmax2]);
%hold on;


plot((tInterp-r/c_inf)*1000+2,pPres*1e-5*r,'Color',cmap(1,:)); %,'LineStyle','--','Color',cmap(i,:));


%% Plot zoom in on peak

handaxes1 = axes('Position',[0.5 0.6 0.35 0.3]);
plot(tdata*1000-tshift,pData*1e-5*r+pshift,'k');

hold on;
xlim([0 10]);

h = plot((output.tPres-r/physConst.c_infty)*1000+2,output.pPres*1e-5*r,...
    'Color',cmap(2,:),'LineStyle','-');
plot((tInterp-r/c_inf)*1000+2,pPres*1e-5*r,'Color',cmap(1,:)); %,'LineStyle','--','Color',cmap(i,:));
