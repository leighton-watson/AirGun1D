%% MAKE FIG 4 SPACE TIME %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Display 1D air gun simulation results. Plot space time plots showing 
% density, temperature, pressure, speed of sound, velocity and Mach
% number inside the air gun chamber

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);

% add code directories
addpath ../SBPSAT
addpath ../sbplib/
addpath ../SeismicAirgunCode

% set color map
colormap(makecmap('dodgerblue',40,20,10));
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
temp = p./(Q.*rho); % temperature 

%% Space-Time Plots %%

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 1200 1000]);

tmin = 0; % minimum time [ms]
tmax = 10; % maximum time [ms]

xmin = 0; % minimum x [m]
xmax = aL; % maximum x [m]

subplot(3,2,1); % density
h = surf(X,T*1000,rho);
view(2); shading interp
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Density');
cb = colorbar; %cb.Label.String = 'kg/m^3';
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent,'XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
h = text(0.01, 9, '(a) density (kg/m^3)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

subplot(3,2,3); % pressure
pa2psi = 0.000145038; % conversion from pa to psi
h = surf(X,T*1000,p*pa2psi);
view(2); shading interp
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Pressure');
cb = colorbar; %cb.Label.String = 'psi';
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent,'XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
h = text(0.01, 9, '(c) pressure (psi)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

subplot(3,2,5); % velocity
h = surf(X,T*1000,v);
view(2); shading interp
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Velocity');
cb = colorbar; %cb.Label.String = 'm/s';
caxis([0 350])
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent,'XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
h = text(0.01, 9, '(e) velocity (m/s)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

p0 = p(1,1); % initial pressure
rho0 = rho(1,1); % initial density
v0 = v(1,1); % initial velocity
c0 = (gamma*p0/rho0)^(0.5); % initial speed of sound
xc = -c0*t + aL; % position of sound wave
zc = ones(size(xc))*1e10; % large number so that sound speed plots on top of space-time figures
xc2 = aL/c0 + c0*(t-aL/c0); % position of reflected sound wave

subplot(3,2,1);
plot3(xc,t*1000,zc,'r--'); 
plot3(xc2,t*1000,zc,'r--');

subplot(3,2,3);
plot3(xc,t*1000,zc,'r--');
plot3(xc2,t*1000,zc,'r--');

subplot(3,2,5);
plot3(xc,t*1000,zc,'r--');
plot3(xc2,t*1000,zc,'r--');

subplot(3,2,2); % temperature
h = surf(X,T*1000,temp);
view(2); shading interp;
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Velocity');
cb = colorbar; caxis([150 300])
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
plot3(xc,t*1000,zc,'r--');
plot3(2*aL/c0-xc,t*1000,zc,'r--');
h = text(0.01, 9, '(b) temperature (K)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

subplot(3,2,4); % speed of sound
h = surf(X,T*1000,c);
view(2); shading interp;
ylabel('Time (ms)'); xlabel('Position (m)'); %title('Velocity');
cb = colorbar; caxis([270 340])
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
plot3(xc,t*1000,zc,'r--'); 
plot3(2*aL/c0-xc,t*1000,zc,'r--');
h = text(0.01, 9, '(d) speed of sound (m/s)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

subplot(3,2,6); % Mach number
h = surf(X,T*1000,v./c);
view(2); shading interp;
ylabel('Time (ms)'); 
xlabel('Position (m)'); %title('Velocity');
cb = colorbar; caxis([0 1])
ylim([tmin tmax]); xlim([xmin xmax]);
set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
hold on;
plot3(xc,t*1000,zc,'r--'); 
plot3(2*aL/c0-xc,t*1000,zc,'r--');
h = text(0.01, 9, '(f) Mach number');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');