%% MAKE FIG 12 Lumped Parameter %%
%
% Make figure for Euler air gun "Geophysics" paper
%
% Compare Euler airgun model with the lumped parameter model

clear all;
clc;
%close all;

addpath ../SBPSAT
addpath ../SeismicAirgunCode

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap_orginal = get(gca,'ColorOrder');
cmap = cmap_orginal; %[cmap_orginal(5,:); cmap_orginal(2,:)];

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 1200]);
alpha = 0.8; % plotting transparency

%%% Parameters %%%

r = 75; % distance from source to receiver [m]
c_inf = 1482; % speed of sound in water [m/s]
rho_inf = 1000; % density in water [kg/m^3]

tmin = 0; % minimum time to plot [ms]
tmax1 = 10; % maximum time to plot [ms]
tmax2 = 500; % second maximum time to plot [ms]

Q = 287.06; % specific gas constant for dry air [J/kgK]
gamma = 1.4; % ratio of heat capacities for dry air
cv = 718; % heat capacity of air at constant volume [J/kgK]
cp = 1010; % heat capacity of air at constant pressure [J/kgK]
T_inf = 288; % ambient water temperature [K]

nx = 100; % number of grid points per 1 m of air gun length

aP = 2000; % air gun pressure [psi]
aA = 12.5; % air gun port area [in^2] % cross-sectional area = port area
aD = 7.5; % air gun depth [m]
m_in = 39.3701; % conversion from m to in
pa_psi = 0.000145038; % conversion from pa to psi
airgunLengths = [1.2 0.6]; % air gun length [m]

%% Iterate over lengths %%

for i = 1:length(airgunLengths)
    
    aL = airgunLengths(i) % display air gun length
    
    %%% Lumped Parameter Model %%%
    
    aV = aL*m_in * aA; % air gun volume [in^3]
    input = [aP, aV, aA]; % inputs for lumped parameter model [pressure, volume, port area]
    physConst = physical_constants(aD,r); % save physical constants. Specific depth and distance from source to receiver
    output = AirgunBubbleSolveOutput(input, physConst, false); % solve lumped parameter model
    
    %%% Plot Lumped Parameter Outputs %%%
    
    % air gun pressure
    subplot(4,1,1);
    h = plot(output.t*1000,output.Y(5,:)*pa_psi,'Color',cmap(i,:),'LineStyle',':');
    h.Color(4) = alpha;
    hold on;
    xlim([tmin tmax1]);
    
    % bubble mass
    subplot(4,1,2);
    h = plot(output.t*1000,output.Y(3,:),'Color',cmap(i,:),'LineStyle',':');
    h.Color(4) = alpha;
    hold on;
    xlim([tmin tmax1]);    
    
    % short time series of pressure perturbation
    subplot(4,1,3);
    h = plot((output.tPres-r/physConst.c_infty)*1000,output.pDir*1e-5*r,...
        'Color',cmap(i,:),'LineStyle',':');
    h.Color(4) = alpha;
    xlim([tmin tmax1]);
    hold on;
    
    % longer time series of pressure perturbation
    subplot(4,1,4);
    h = plot((output.tPres-r/physConst.c_infty)*1000,output.pDir*1e-5*r,...
        'Color',cmap(i,:),'LineStyle',':');
    h.Color(4) = alpha;
    xlim([tmin tmax2]);
    hold on;
    
    
    
    
     
    %%% Euler Air Gun Model %%%
    
    sol = runEulerCode_initBubbleVol(nx, aV, aP, aL, aA, aD);
    
    t = sol.x; % time
    x = [0:ceil(aL*nx)]./nx; % space vector
    [T,X] = meshgrid(t,x); % create mesh for space-time plots
    
    % initialize matrices
    rho = zeros(length(x), length(t));
    rhov = zeros(length(x), length(t));
    e = zeros(length(x), length(t));
    
    for j = 1:length(x) % extract air gun properties
        rho(j,:) = sol.y(3*j+2,:); % density
        rhov(j,:) = sol.y(3*j+3,:); % density * velocity
        e(j,:) = sol.y(3*j+4,:); % internal energy
    end
    
    v = rhov./rho; % velocity [m/s]
    p = (gamma-1)*(e-0.5*rho.*v.^2); % pressure
    
    % bubble properties
    R = sol.y(1,:); % bubble radius [m]
    U = sol.y(2,:); % bubble wall velocity [m/s]
    V = 4/3*pi.*R.^3; % bubble volume [m^3]
    m = sol.y(3,:); % bubble mass [kg]
    E = sol.y(4,:); % internal energy [J]
    Temp = E./(m.*cv); % temperature [K]
    pb = m.*Q.*Temp./V; % bubble pressure [Pa]
    [~, solDY] = deval(sol, t);
    A = solDY(2,:); % acceleration
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    
    %%% Plot Euler Air Gun Outputs %%%
    
    % air gun pressure
    subplot(4,1,1);
    idx = length(x)-1; % spatial position to plot
    plot(t*1000,p(idx,:)*pa_psi,'Color',cmap(i,:),'LineStyle','-');
    
    
    % bubble mass
    subplot(4,1,2);
    plot(t*1000,m,'Color',cmap(i,:),'LineStyle','-');
    
    % short time series of pressure perturbation
    subplot(4,1,3);
    plot((tDir-r/c_inf)*1000, pDir*1e-5*r,'Color',cmap(i,:),'LineStyle','-');
    
    % longer time series of pressure perturbation
    subplot(4,1,4);
    plot((tDir-r/c_inf)*1000, pDir*1e-5*r,'Color',cmap(i,:),'LineStyle','-');    
    

end

%% Format Figures %%

figure(1);
subplot(4,1,1);
ylim([0 2400]);
ylabel('psi');
h = text(0.2,2200,'(a) air gun pressure');
set(h,'FontSize',24);
set(h,'FontWeight','bold');

subplot(4,1,2);
ylabel('kg');
h = text(0.2,1.35,'(b) bubble mass');
set(h,'FontSize',24);
set(h,'FontWeight','bold');

subplot(4,1,3);
ylabel('bar m');
ylim([0 8]);
h = text(0.2,7.2,'(c) pressure perturbation');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
xlabel('Time (ms)');

subplot(4,1,4);
ylabel('bar m');
ylim([-2 8]);
h = text(10,7,'(d) pressure perturbation');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
xlabel('Time (ms)');
legend('L = 2 m','L = 0.5 m');
 