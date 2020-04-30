%% MAKE FIG A2 NUMERICAL RESOLUTION %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Run simulations for different grid resolution to demonstrate numerical
% resolution of 1D air gun simulations.

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

% add directories
addpath ../sbplib/
addpath ../SBPSAT
addpath ../SeismicAirgunCode

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 700]);

%% 1D Air Gun Code %%

gamma = 1.4; % ratio of heat capacities
cv = 718; % heat capacity of air at constant volume [J/kgK]
cp = 1010; % heat capacity of air at constant pressure [J/kgK]
Q = 287.06; % specific gas constant for dry air [J/kgK]
    
r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

pa2psi = 0.000145038; % conversion from pa to psi

aP = 2000; % air gun pressure [psi]
aL = 1; % air gun length [m]
aA = 16; % air gun port area [in^2] % cross-sectional area = port area
aD = 7.5; % air gun depth [m]

nx = [10 100 1000]; % number of grid points per 1 m of air gun length
nx_plot = [10 100 1000]; % grid points to plot outlet pressure

Rsave = []; % save bubble radius
Psave = []; % save acoustic pressure

for j = 1:length(nx)
    
    nx(j)
    sol = runEulerCode(nx(j), aP, aL, aA, aD);

    %%% Save Outputs %%%
    
    t = sol.x; % time
    x = [0:ceil(aL*nx(j))]./nx(j); % space vector
    idx = length(x)-1; % spatial position to plot
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
    
    v = rhov./rho; % velocity [m/s]
    p = (gamma-1)*(e-0.5*rho.*v.^2); % pressure
    c = (gamma*p./rho).^(0.5); % speed of sound
    
    % outlet pressure
    if ismember(nx(j), nx_plot)
        figure(1);
        subplot(2,1,1);
        plot(t*1000, p(idx,:)*pa2psi);
    hold on;
    
    end
    
    % bubble properties
    R = sol.y(1,:); % bubble radius [m]
    U = sol.y(2,:); % bubble wall velocity [m/s]
    m = sol.y(3,:); % bubble mass [kg]
    [~,solDY] = deval(sol, t);
    A = solDY(2,:); % acceleration
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    
    % compute peak pressure, slope and rise time
    [ppeak(j),idx2] = max(pDir*1e-5*r);
    slope(j) = ppeak(j)/((tDir(idx2)-r/c_inf)*1000);
    riseTime(j) = ((tDir(idx2)-r/c_inf)*1000);
       
end

%% Plotting %%

figure(1);
subplot(2,1,1);
p_analytical = aP * (2/(gamma+1))^((2*gamma)/(gamma-1));
plot([min(t) max(t)],[p_analytical p_analytical],'k--');
ylim([0 aP]);
xlim([0 1]);
xlabel('Time (ms)');
ylabel('psi');

figure(1);
subplot(2,1,2);
semilogx(nx, ppeak,'ko:','markerfacecolor','k');
hold on;
plot(nx, slope,'ko--','markerfacecolor','k');
plot(nx, riseTime,'ko-','markerfacecolor','k');
xlabel('Number of grid points');
xlim([min(nx) max(nx)]);
ylim([0 7]);
grid on;
legend('peak acoustic pressure (bar m)','slope (bar m/ms)','rise time (ms)');


