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

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 260]);

figure(1); clf;
figHand2 = figure(2); clf;
set(figHand2,'Position',[100 100 600 900]);


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

%aP = [100 220 420 610 810 1030 1200 1400 1600 1800 2000]; % air gun pressure for slope/rise time plot
%aP = [1000 2000]%
aP_plot = [220 420 610 810 1030]; % air gun pressures to plot [psi]
%aP_plot = aP;
%aP_plot = [220 420 610 810 1030];
%aP_plot = [220 1030];
%aP_plot = [2000 1000]
aP_plot = fliplr(aP_plot);
aP = aP_plot;

%aL = 0.6; % air gun length [m]
%aA = 16; % air gun port area [in^2] % cross-sectional area = port area
aL = 1.2;
aA = 12.5;
aD = 7.5; % air gun depth [m]
aV = aL*39.3701*aA;

physConst = physical_constants(aD,r);
k = [5 4 3 2 1];
j = 1;
for i = 1:length(aP)
    
    aP(i)
    
    sol = runEulerCode(nx, aP(i), aL, aA, aD);
    
    t = sol.x; % time
    
    % bubble properties
    R = sol.y(1,:); % bubble radius [m]
    U = sol.y(2,:); % bubble wall velocity [m/s]
    m = sol.y(3,:); % bubble mass [kg]
    [~,solDY] = deval(sol, t);
    A = solDY(2,:); % acceleration
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    
    if ismember(aP(i), aP_plot)
        
        % bubble mass
        figure(2);
        subplot(3,1,1);
        h = plot(t*1000, m,'Color',cmap(j,:));
        %h.Color(4) = 0.6;
        hold on;
        xlim([tmin tmax]);
        ylabel('kg');
        %ylim([0 0.6]);
        %xlabel('Time (ms)');
        
        
        % plot analytical mass flow rate
        dmdt = (aP(i)*psi_pa)*(aA*in2_m2) * (gamma/(Q*T_inf))^(1/2) * ...
            (2/(gamma+1))^((gamma+1)/(gamma-1));
        m_analytical = m(1) + t*dmdt;
        h = plot(t*1000,m_analytical,'Color',cmap(j,:),'LineStyle','--');
        h.Color(4) = 0.8;
        
        
        % pressure perturbation
        subplot(3,1,2);
        plot((tDir-r/c_inf)*1000,pDir*1e-5*r,'Color',cmap(j,:));
        hold on;
        xlim([tmin tmax]);
        ylabel('bar m');
        xlabel('Time (ms)');
        
        % analytical acoustic pressure perturbation
        input = [aP(i), aV, aA];
        airgunInit = airgun_initialization(input, physConst);
        bubbleInit = bubble_initialization(airgunInit, physConst);
        initCond = initCond_save(bubbleInit, airgunInit);
        initCond(7) = dmdt;
        params = params_save(bubbleInit, airgunInit, physConst);
        options = odeset('RelTol',1e-6);
        sol = ode45(@modified_herring_eqn_constmass, [0 10], initCond, options, params);
        t2 = sol.x'; %time vector
        R2 = sol.y(1,:)'; %bubble radius
        U2 = sol.y(2,:)'; %bubble wall velocity
        
        m2 = sol.y(3,:)'; % bubble mass
        subplot(3,1,1);
        plot(t2*1000, m2,'Color',cmap(j,:),'LineStyle',':');
        
        subplot(3,1,2);        
        [~, solDY] = deval(sol, t2);
        A2 = solDY(2,:)'; %bubble wall acceleration
       
        [tDir2, pDir2] = pressure_eqn(t2, R2, U2, A2, physConst.rho_infty, physConst.c_infty, physConst.r); %direct arrival
        plot((tDir2-r/c_inf)*1000-0.7933, pDir2*1e-5*r,'Color',cmap(j,:),'LineStyle','--');
        
        %%% longer time series including ghost
        figure(1);
        [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, ... %ghost
        r + 2*aD);
    dt = 1e-5;
    tInterp = min(tDir):dt:max(tDir);
    pDirInterp = pchip(tDir, pDir, tInterp);
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp;
        
        plot((tInterp-r/c_inf)*1000,pPres*1e-5*r,'Color',cmap(k(j),:));
        hold on;
        xlim([tmin 300]);
        ylim([-1.5 3])
        ylabel('bar m');
        xlabel('Time (ms)');
        
        j = j+1;
    end
    
    % compute peak pressure, slope and rise time
    [ppeak(i),idx] = max(pDir*1e-5*r);
    slope(i) = ppeak(i)/((tDir(idx)-r/c_inf)*1000);
    riseTime(i) = ((tDir(idx)-r/c_inf)*1000);
 
end

%% 
figure(2);
subplot(3,1,3);
plot(aP, riseTime,'k-');
hold on;
plot(aP, slope,'k--');
plot(aP, ppeak,'k:');
xlim([min(aP), max(aP)]);
xlabel('Air gun pressure (psi)');
legend('Rise Time (ms)','Slope (bar m/ms)','Peak Pressure (bar m)');
ylim([0 7])
h = text(150,5.2,'(c)');
set(h,'FontSize',24);
set(h,'FontWeight','bold');

subplot(3,1,1);
ylim([0 0.8]);
h = text(0.2,0.52,'(a) bubble mass');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
    
subplot(3,1,2);
ylim([0 3])
h = text(0.2,3.5,'(b) acoustic pressure');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
%legend('220 psi','420 psi','610 psi','810 psi','1030 psi')