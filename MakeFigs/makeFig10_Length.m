%% LENGTH %%


clear all;
clc;
%close all;

addpath ../SBPSAT
addpath ../SeismicAirgunCode

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

figHand2 = figure(2); clf;
set(figHand2,'Position',[100 100 600 600]);


%% Run Euler Air Gun Simulation %%

nx = 50; % number of grid points per 1 m of air gun length

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


aL_plot = [2 1 0.6 0.4 0.2];
%aL = aL_plot;
aL = [0.2 0.3 0.4 0.5 0.6 0.8 1 1.2 1.4 1.6 1.8 2];
aL = fliplr(aL);

aP = 1020;
aA = 12.5;
aD = 7.5;
j = 1;
for i = 1:length(aL)
    
    aL(i)
    
    sol = runEulerCode(nx, aP, aL(i), aA, aD);
    
    t = sol.x; % time
    
    % bubble properties
    R = sol.y(1,:); % bubble radius [m]
    U = sol.y(2,:); % bubble wall velocity [m/s]
    m = sol.y(3,:); % bubble mass [kg]
    [~,solDY] = deval(sol, t);
    A = solDY(2,:); % acceleration
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    
    if ismember(aL(i), aL_plot)
        
        % pressure perturbation
        subplot(2,1,1);
        plot((tDir-r/c_inf)*1000,pDir*1e-5*r);
        hold on;
        xlim([tmin tmax]);
        ylabel('bar m');
        xlabel('Time (ms)');
        
        
    end
    
    % compute peak pressure, slope and rise time
    [ppeak(i),idx] = max(pDir*1e-5*r);
    slope(i) = ppeak(i)/((tDir(idx)-r/c_inf)*1000);
    riseTime(i) = ((tDir(idx)-r/c_inf)*1000);
 
end

%% 

subplot(2,1,2);
plot(aL, riseTime,'k-');
hold on;
plot(aL, slope,'k--');
plot(aL, ppeak,'k:');
xlim([min(aL), max(aL)]);
xlabel('Air Gun Length (m)');
legend('Rise time (ms)','Slope (bar m/ms)','Peak pressure (bar m)');
ylim([0 5])
% h = text(150,5.2,'(c)');
% set(h,'FontSize',24);
% set(h,'FontWeight','bold');

subplot(2,1,1);
h = text(0.2,2.6,'(a) acoustic pressure');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
legend('2 m','1 m','0.6 m','0.4 m','0.2 m')
ylim([0 3])