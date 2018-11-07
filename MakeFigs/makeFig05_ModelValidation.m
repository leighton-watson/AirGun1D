%% MAKE FIG 4 Model Validation %%
%
% Make figure for Euler air gun "Geophysics" paper
%
% Simulate air gun dynamics. Plot time snapshots of the internal dynamics.
% Compare numerical and analytical solution

clear all;
clc;
%close all;

addpath ../SBPSAT
addpath ../SeismicAirgunCode

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
colormap(makecmap('black',40,10,10));
%colormap([makecmap('orangered',50,20,20);flipud(makecmap('tomato',50,20,20))]);
cmap = get(gca,'ColorOrder');

%% Run Euler Air Gun Simulation %%

nx = 50; % number of grid points per 1 m of air gun length

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

%% Time Snapshots %%

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 1200]);

p0 = p(1,1); % initial pressure
rho0 = rho(1,1); % initial density
v0 = v(1,1); % initial velocity
c0 = (gamma*p0/rho0)^(0.5); % initial speed of sound

pa2psi = 0.000145038; % conversion from pa to psi

v_analytical = zeros(size(x));
c_analytical = zeros(size(x));
p_analytical = zeros(size(x));
rho_analytical = zeros(size(x));

alpha = 0.4; % plotting transparency

tplot = [1 2 3 4]; % times to plot [ms]

for i = 1:length(tplot)
    
    [~,idx] = min(abs(t-tplot(i)/1000)); % find nearest value
    
    %%% analytical solutions %%%
    position_expansionfan = -c0*t(idx); % position of front of expansion fan
    for j = 1:length(x)
        if x(j) < position_expansionfan+aL % undisturbed properties
            v_analytical(j) = v0;
            c_analytical(j) = c0;
            p_analytical(j) = p0;
            rho_analytical(j) = rho0;
        else % properties inside expansion fan
            % v_analytical(j) = (2/(gamma+1)).*((x(j)-aL)./t(idx)+((gamma-1)/2)*v0+c0);
            % c_analytical(j) = v_analytical(j)-((x(j)-aL))/t(idx);
            
            %v_analytical(j) = (2/(gamma+1)).*((aL-x(j))./t(idx)+((gamma-1)/2)*v0+c0);
            v_analytical(j) = (2/(gamma+1)).*(c0-(aL-x(j))/t(idx));
            c_analytical(j) = v_analytical(j)+((aL-x(j)))/t(idx);
            
            p_analytical(j) = p0*(c_analytical(j)/c0)^((2*gamma)/(gamma-1));
            rho_analytical(j) = (gamma*p_analytical(j))/c_analytical(j)^2;
        end
    end
    
    subplot(4,1,1); % density
    figure(1);
    h = plot(x,rho_analytical,'-','Color',cmap(i,:));
    h.Color(4) = alpha;
    ylabel('kg/m^3'); %title('Density');
    hold on; grid on;
    ylim([50 220])
    xlim([0 aL]);
    set(h.Parent,'YTick',[50 100 150]);
    set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
    h = text(0.01, 200,'(a) density');
    set(h,'FontSize',24);
    set(h,'FontWeight','bold');
    
    % h = vline(1.2-t(idx)*c0);
    h.Color = cmap(i,:);
    
    
    subplot(4,1,2); % pressure
    h = plot(x,p_analytical*pa2psi,'-','Color',cmap(i,:));
    ylabel('psi'); %title('Pressure');
    h.Color(4) = alpha;
    hold on; grid on;
    ylim([500 2500])
    xlim([0 aL]);
    set(h.Parent,'YTick',[500 1000 1500 2000]);
    set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
    h = text(0.01, 2300,'(b) pressure');
    set(h,'FontSize',24);
    set(h,'FontWeight','bold');
    
    
    subplot(4,1,3); % velocity
    h = plot(x,v_analytical,'-','Color',cmap(i,:));
    ylabel('m/s'); %title('Velocity');
    h.Color(4) = alpha;
    hold on; grid on;
    ylim([-30 300])
    xlim([0 aL]);
    set(h.Parent,'YTick',[0 100 200 300]);
    set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
    h = text(0.01, 265,'(c) velocity');
    set(h,'FontSize',24);
    set(h,'FontWeight','bold');
    
    
    subplot(4,1,4); % speed of sound
    h = plot(x,c_analytical,'-','Color',cmap(i,:));
    h.Color(4) = alpha;
    ylabel('m/s'); xlabel('Position (m)'); %title('Speed of Sound');
    hold on; grid on;
    ylim([260 370])
    xlim([0 aL]);
    set(h.Parent,'YTick',[260 300 340]);
    set(h.Parent','XTick',[0 0.2 0.4 0.6 0.8 1 1.2]);
    h = text(0.01, 358,'(d) speed of sound');
    set(h,'FontSize',24);
    set(h,'FontWeight','bold');
    
end

for i = 1:length(tplot)
    
    [~,idx] = min(abs(t-tplot(i)/1000)); % find nearest value
    
    %%% numerical solutions %%%
    subplot(4,1,1); % density
    h = plot(x,rho(:,idx),'x','Color',cmap(i,:));

    subplot(4,1,2); % pressure
    plot(x,p(:,idx)*pa2psi,'x','Color',cmap(i,:));

    subplot(4,1,3); % velocity 
    plot(x,v(:,idx),'x','Color',cmap(i,:));
    
    subplot(4,1,4); % speed of sound
    h = plot(x,c(:,idx),'x','Color',cmap(i,:));
       
end
% 
% subplot(4,1,4);
% legend('t = 1 ms','t = 2 ms','t = 3 ms','t = 4 ms','Location','SouthWest')
% 
% print -depsc 'Fig4_ModelValidation';