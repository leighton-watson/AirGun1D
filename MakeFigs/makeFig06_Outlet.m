%% MAKE FIG 6 OUTLET %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Display 1D air gun simulation results. Plot properties at outlet and
% compare to analytical solutions.

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);

% add code directories
addpath ../SBPSAT
addpath ../SeismicAirgunCode

cmap = get(gca,'ColorOrder');

%% Run Euler Air Gun Simulation %%

nx = 100; % number of grid points per 1 m of air gun length

aP = 2000; % air gun pressure [psi]
aL = 1.2; % air gun length [m]
aA = 12.5; % air gun port area [in^2] % cross-sectional area = port area
aD = 7.5; % air gun depth [m]

sol = runEulerCode(nx, aP, aL, aA, aD); 

%% Save Outputs %%

t = sol.x; % time
x = [0:ceil(aL*nx)]./nx; % space vector
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
cv = 718; % heat capacity of air at constant volume [J/kgK]
cp = 1010; % heat capacity of air at constant pressure [J/kgK]
Q = 287.06; % specific gas constant for dry air [J/kgK]

v = rhov./rho; % velocity [m/s]
p = (gamma-1)*(e-0.5*rho.*v.^2); % pressure
c = (gamma*p./rho).^(0.5); % speed of sound

% bubble properties
R = sol.y(1,:); % bubble radius [m]
U = sol.y(2,:); % bubble wall velocity [m/s]
V = 4/3*pi.*R.^3; % bubble volume [m^3]
m = sol.y(3,:); % bubble mass [kg]
E = sol.y(4,:); % internal energy [J]
Temp = E./(m.*cv); % temperature [K]
pb = m.*Q.*Temp./V; % bubble pressure [Pa]

p0 = p(1,1); % initial pressure
rho0 = rho(1,1); % initial density
v0 = v(1,1); % initial velocity
c0 = (gamma*p0/rho0)^(0.5); % initial speed of sound

pa2psi = 0.000145038; % conversion from pa to psi

%% Properties at Outlet %%

% compute analytical expressions for properties at the outlet
ve = 2/(gamma+1)*c0; % exit velocity
ce = ve; % exit sound speed
pe = p0*(2/(gamma+1))^((2*gamma)/(gamma-1)); % exit pressure
rhoe = (gamma*p0*(2/(gamma+1))^((2*gamma)/(gamma-1)))/((2/(gamma+1)*c0)^2); % exit density


figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 900]);
tmin = 0; tmax = 10;
idx = length(x)-1; % spatial position to plot

subplot(3,1,1); % density
h = plot(t*1000,rho(idx,:),'Color',cmap(1,:));
xlim([tmin tmax]);
ylabel('kg/m^3'); %xlabel('Time (s)'); %title('Density');
hold on; %grid on;
plot([tmin tmax],[rhoe rhoe],'k--');
h = text(0.2, 180, '(a) density');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
h2 = vline(2*aL/c0*1000);

subplot(3,1,2); % pressure
plot(t*1000,p(idx,:)*pa2psi,'Color',cmap(1,:));
xlim([tmin tmax]); ylim([0 2000])
ylabel('psi'); %xlabel('Time (s)'); %title('Pressure');
hold on; %grid on;
h = text(0.2, 1800, '(b) pressure');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
h2 = vline(2*aL/c0*1000);

% ambient pressure
rho_inf = 1e3; % density [kg/m^3]
pa = 1e5; % atmospheric pressure [Pa]
g = 9.8; % gravitational acceleration [m/s^2]
p_inf = pa + rho_inf*g*aD; % ambient pressure at depth [Pa]
plot([min(t) max(t)]*1000,[p_inf,p_inf]*pa2psi,'Color',cmap(3,:),'LineStyle','--')

% bubble pressure
plot(t*1000, pb*pa2psi,'Color',cmap(2,:),'LineStyle',':');
legend('Air gun pressure','Ambient pressure','Bubble pressure');
plot([tmin tmax],[pe*pa2psi pe*pa2psi],'k--');

subplot(3,1,3); % velocity
plot(t*1000,v(idx,:),'Color',cmap(1,:));
xlim([tmin tmax]);
ylabel('m/s'); %title('Velocity');
xlabel('Time (s)');
hold on; %grid on;
ylim([0 400]);
h = text(0.2, 360, '(c) velocity');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
h2 = vline(2*aL/c0*1000);

% speed of sound
plot(t*1000, c(idx,:),'--','Color',cmap(6,:));
legend('Air gun velocity','Speed of sound','Location','SouthEast');
plot([tmin tmax],[ve ve],'k--');