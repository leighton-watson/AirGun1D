%% DATA %%

clear all;

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 1200]);

cmap = get(gca,'ColorOrder');
msize = 8;

%% TIME SERIES %%%
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/598ci/FarField'
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/50ci/FarField'

subplot(4,1,1);
hold on;

m = 32000;
dt = 31.25e-6;
t = 0:dt:(m-1)*dt;
t = t-40/1000;

tshift = [51.28 51.72 50.56 50.69 51.97];
pshift = [0.003751 0.0006554 -0.002296 -0.006076 0.001113];

% 1030 psi
data = load('188_0750cm_1030psi_598ci_DHA.csv');
plot(t*1000-tshift(1), data(:,1)*(1e-5*75)-pshift(1));

% 810 psi
data = load('183_0750cm_0810psi_598ci_DHA.csv');
plot(t*1000-tshift(2), data(:,1)*(1e-5*75)-pshift(2));

% 610 psi
data = load('174_0750cm_0610psi_598ci_DHA.csv');
plot(t*1000-tshift(3), data(:,1)*(1e-5*75)-pshift(3));

% 420 psi
data = load('173_0750cm_0420psi_598ci_DHA.csv');
plot(t*1000-tshift(4), data(:,1)*(1e-5*75)-pshift(4));

% 220 psi
data = load('166_0750cm_0220psi_598ci_DHA.csv');
plot(t*1000-tshift(5), data(:,1)*(1e-5*75)-pshift(5));

xlim([-10 200]);
ylim([-1 1.5]);
xlabel('Time (ms)');
ylabel('bar m');
box on


%% RISE TIME %%

subplot(4,1,2);
hold on;

addpath '/Users/lwat054/Documents/Stanford_University/Papers/Submitted/Low Frequency Pneumatic Seismic Sources/Figs/RiseTime_Period/RiseTime'

% 25 m depth
data = load('rtim_2500cm_50ci_P-rtim.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(5,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));

data = load('rtim_2500cm_598ci_P-rtim.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(5,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));

% 15 m depth
% data = load('rtim_1500cm_50ci_P-rtim.txt');
% plot(data(:,1), data(:,2),'^','Color',cmap(4,:),'MarkerSize',msize,...
%     'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor',cmap(4,:));

data = load('rtim_1500cm_598ci_P-rtim.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(4,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor',cmap(4,:));

% 10 m depth
data = load('rtim_1000cm_50ci_P-rtim.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(3,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));

data = load('rtim_1000cm_598ci_P-rtim.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(3,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));

% 7.5 m depth
data = load('rtim_0750cm_50ci_P-rtim.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(2,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));

data = load('rtim_0750cm_598ci_P-rtim.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(2,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));

% 5m depth
data = load('rtim_0500cm_50ci_P-rtim.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(1,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));

data = load('rtim_0500cm_598ci_P-rtim.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(1,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));



% label
xlabel('Pressure (psi)');
ylabel('Rise Time (ms)');
box on;
xlim([200 1200]);
ylim([1 5.5])



%% PEAK AMPLITUDE %%

subplot(4,1,3);
hold on;


%%% 598 ci %%%
load data_attributes.mat

for i = 1:length(csvFiles_598)
    if depth(i) == 2500
        plot(pressure(i), ppeak(i),'o','Color',cmap(5,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));
        
    elseif depth(i) == 1500
        plot(pressure(i), ppeak(i),'o','Color',cmap(4,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor',cmap(4,:));
        
    elseif depth(i) == 1000
        plot(pressure(i), ppeak(i),'o','Color',cmap(3,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));
        
    elseif depth(i) == 750
        plot(pressure(i), ppeak(i),'o','Color',cmap(2,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));
        
    elseif depth(i) == 500
        plot(pressure(i), ppeak(i),'o','Color',cmap(1,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));
    end
end


for i = 1:length(csvFiles_50)
    if depth2(i) == 2500
        plot(pressure2(i), ppeak2(i),'^','Color',cmap(5,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));
        
    elseif depth2(i) == 1500
        plot(pressure2(i), ppeak2(i),'^','Color',cmap(4,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor',cmap(4,:));
        
    elseif depth2(i) == 1000
        plot(pressure2(i), ppeak2(i),'^','Color',cmap(3,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));
        
    elseif depth2(i) == 750
        plot(pressure2(i), ppeak2(i),'^','Color',cmap(2,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));

    elseif depth2(i) == 500
        plot(pressure2(i), ppeak2(i),'^','Color',cmap(1,:),'MarkerSize',msize,...
            'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));
    end
end

% label
xlabel('Pressure (psi)');
ylabel('Peak Amplitude (bar m)');
box on;
xlim([200 1200]);
ylim([0 1]);

%% SLOPE %%

subplot(4,1,4);
hold on;

addpath '/Users/lwat054/Documents/Stanford_University/Papers/Submitted/Low Frequency Pneumatic Seismic Sources/Figs/Pressure/PressureSlopeData'

% 25 m depth
data = load('rtim2_2500cm_50ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(5,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));

data = load('rtim2_2500cm_598ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(5,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor',cmap(5,:));

% 15 m depth
data = load('rtim2_1500cm_598ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(4,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor',cmap(4,:));

% 10 m depth
data = load('rtim2_1000cm_50ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(3,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));

data = load('rtim2_1000cm_598ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(3,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor',cmap(3,:));

% 7.5 m depth
data = load('rtim2_0750cm_50ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(2,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));

data = load('rtim2_0750cm_598ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(2,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor',cmap(2,:));

% 5m depth
data = load('rtim2_0500cm_50ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'^','Color',cmap(1,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));

data = load('rtim2_0500cm_598ci_Pressure-Slope.txt');
plot(data(:,1), data(:,2),'o','Color',cmap(1,:),'MarkerSize',msize,...
    'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor',cmap(1,:));


% label
xlabel('Pressure (psi)');
ylabel('Slope (bar m/ms)');
xlim([200 1200])
ylim([0 0.5])
box on;