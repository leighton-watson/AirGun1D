%% COMPARE BUBBLE MODELS %%
%
% Compare between the Herring, Modified Herring, and Gilmore equation.
% Should I include the Rayleigh equation???
%
% Generate the figure for the appendix of the seismic airgun paper.

clear all; clc; 
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',22);
cmap = get(gca,'ColorOrder');

%% INITIAL PARAMETERS %%

%%%TIME%%%
time = [0 10]; %time interval to intergrate over [s]
time_plot = [0 0.5]; % time interval to plot over [s]

%%%INITIAL CONDITIONS%%%
init_pressure_psi = 2000; 
init_pressure = init_pressure_psi * 6894.8;

init_volume_in3 = 400;
init_volume = init_volume_in3 * 1.63871e-5;
init_radius = ((3*init_volume)/(4*pi))^(1/3);
init_vel = 0;

T_inf = 288; % ambient temperature
Q = 287.06; % gas constant
init_mass = init_pressure*init_volume/(T_inf*Q); % mass in bubble

init_cond = [init_radius, init_vel, T_inf, init_pressure]; %initial conditions for solver


%%%PARAMETERS%%%
rho_infty = 1000; %ambient density [kg/m^3]
p0 = 1e5; %atmospheric pressure [Pa]
depth = 7; % depth of air gun [m]
g = 9.8; % gravitational acceleration [m/s^2]
p_infty = p0 + rho_infty*g*depth; % pressure at depth of air gun
k = 1.4; %polytropic index, for isothermal k = 1, for adiabatic k = 1.4
B = 3.04913e8; %constant from Tait equation of state [Pa]
n = 7.15; %constant from Tai equation of state [ ] 
c_infty = sqrt(n*(p_infty+B)/rho_infty); %speed of sound in fluid [m/s]

params = [rho_infty, p_infty, init_pressure, init_radius, k, c_infty, init_mass, B, n]; %parameters for solver

%%%SOLVER OPTIONS%%%
options = odeset('RelTol',1e-4); %set options for ODE solver

%%%FIGURE%%%
figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 600]);

%% SOLVER %%

%%MODIFIED HERRING%%%
[tmh,Ymh] = ode45(@modifiedHerring, time, init_cond, options, params);
subplot(2,1,1);
h2 = plot(tmh*1000, Ymh(:,1),'Color',cmap(1,:));
ylabel('m');
%grid on;
hold on;
plot([time(1) time(end)]*1000,[Ymh(end,1) Ymh(end,1)],'Color',cmap(1,:),'LineStyle','--');
xlim([time_plot(1) time_plot(end)]*1000);
ylim([0 1]);
h = text(10,0.9,'(a) Bubble radius');
set(h,'FontSize',24);
set(h,'FontWeight','bold');



%%%HERRING%%%
[th,Yh] = ode45(@Herring, time, init_cond, options, params);
%plot(th*1000, Yh(:,1));

%%%GILMORE%%%
[tg,Yg] = ode45(@Gilmore, time, init_cond, options, params);
%plot(tg*1000, Yg(:,1));

%%%RAYLEIGH%%%
[tR,YR] = ode45(@Rayleigh, time, init_cond, options, params);
%plot(tR*1000, YR(:,1));
%% DIFFERENCE BETWEEN BUBBLE MODELS %%

% interpolate on to a common time vector
Yh_int = pchip(th,Yh(:,1),tmh);
Yg_int = pchip(tg,Yg(:,1),tmh);
YR_int = pchip(tR,YR(:,1),tmh);

% difference between Herring and Modified Herring equation
diff_h = (Yh_int-Ymh(:,1))./(Ymh(:,1)).*100;

% difference between Gilmore and Modified Herring equation
diff_g = (Yg_int-Ymh(:,1))./(Ymh(:,1)).*100;

% difference between Rayleigh and Modified Herring equation
%diff_R = (YR_int-Ymh(:,1))./(Ymh(:,1)).*100;

subplot(2,1,2);
plot(tmh*1000, diff_h,'Color',cmap(2,:)); hold on;
plot(tmh*1000, diff_g,'Color',cmap(3,:));
%plot(tmh*1000, diff_R,'Color',cmap(4,:));
xlim([time_plot(1) time_plot(end)]*1000);
ylim([-0.8 1.8]);
%grid on;
leg = legend('Herring','Gilmore','Location','NorthEast');
set(leg,'FontSize',18);
xlabel('Time (ms)');
ylabel('%')
h = text(10,1.55,'(b) Difference');
set(h,'FontSize',24);
set(h,'FontWeight','bold');

% 
% figHand2 = figure(2); clf;
% set(figHand2,'Position',[200 250 600 450]);
% 
% subplot(2,1,1);
% h2 = plot(tmh, Ymh(:,1),'Color',cmap(2,:));
% title('Modified Herring');
% ylabel('R (m)');
% grid on;
% xlim([time_plot(1) time_plot(end)]);
% 
% % interpolate on to a common time vector
% Yh_int = pchip(th,Yh(:,1),tmh);
% Yg_int = pchip(tg,Yg(:,1),tmh);
% 
% % difference between Herring and Modified Herring equation
% diff_h = (Yh_int-Ymh(:,1)); %./(Ymh(:,1)).*100;
% 
% % difference between Gilmore and Modified Herring equation
% diff_g = (Yg_int-Ymh(:,1)); %./(Ymh(:,1)).*100;
% 
% subplot(2,1,2);
% plot(tmh, diff_h,'Color',cmap(1,:));
% hold on;
% plot(tmh, diff_g,'Color',cmap(3,:));
% xlim([time_plot(1) time_plot(end)])

