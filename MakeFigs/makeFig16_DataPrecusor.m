%% DATA %%

clear all;

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);
cmap = get(gca,'ColorOrder');

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 600 300]);

cmap = get(gca,'ColorOrder');
msize = 8;

%% TIME SERIES %%%
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/598ci/FarField'
addpath '/Users/lwat054/Documents/Stanford_University/Research/SeismicAirguns/Data/Lake/CSVFormat/50ci/FarField'

hold on;

m = 32000;
dt = 31.25e-6;
t = 0:dt:(m-1)*dt;
t = t-40/1000;

tshift = [51.28 51.72 50.56 50.69 51.97];
pshift = [0.003751 0.0006554 -0.002296 -0.006076 0.001113];

% 1030 psi
data = load('188_0750cm_1030psi_598ci_DHA.csv');
plot(t*1000-tshift(1), data(:,1)*(1e-5*75)-pshift(1),'Color',cmap(5,:));

% 810 psi
data = load('183_0750cm_0810psi_598ci_DHA.csv');
plot(t*1000-tshift(2), data(:,1)*(1e-5*75)-pshift(2),'Color',cmap(4,:));

% 610 psi
data = load('174_0750cm_0610psi_598ci_DHA.csv');
plot(t*1000-tshift(3), data(:,1)*(1e-5*75)-pshift(3),'Color',cmap(3,:));

% 420 psi
data = load('173_0750cm_0420psi_598ci_DHA.csv');
plot(t*1000-tshift(4), data(:,1)*(1e-5*75)-pshift(4),'Color',cmap(2,:));

% 220 psi
data = load('166_0750cm_0220psi_598ci_DHA.csv');
plot(t*1000-tshift(5), data(:,1)*(1e-5*75)-pshift(5),'Color',cmap(1,:));

xlim([-3 12]);
ylim([0 1]);
xlabel('Time (ms)');
ylabel('bar m');
box on


