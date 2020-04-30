%% MAKE FIG 1 DATA CWT %%
%
% Watson, Werpers and Dunham (2018) What controls the initial peak of an
% air gun source signature, Geophysics
%
% Display data. Plot acoustic pressure time series and continuous wavelet 
% transform of field data from Lake Seneca.
%
% For information about the data see Ronen and Chelminski (2018) A next 
% generation seismic source with low frequency signal and low 
% environmental impact, 80th EAGE Conference & Exhibition. 
% doi:10.3997/2214-4609.201800745

clear all; clc;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',24);

% add data directory
addpath ../Data
addpath ../SBPSAT/

% specify data
dataStr = '188_0750cm_1030psi_598ci_DHA.csv';

% format figure for plotting
figHand1 = figure(1); clf; 
set(figHand1,'Position',[100 100 600 700]);

% load data
data = csvread(dataStr);
pData = data(:,1);
k = length(data);
dt = 31.25e-6;
Fs = 1/dt;
t = 0:dt:(k-1)*dt;
r = 75;
tshift = -41.25;
pDataBarM = pData*1e-5*r;

%%% plot time series %%%
subplot(3,1,1);
h = plot((t)*1000+tshift,pDataBarM);
set(h.Parent,'XTick',[0 100 200 300 400 500]);
xlim([0 500])
ylabel('\Delta p (bar m)');
h = text(0.006*1000, 0.74, '(a)');
set(h,'FontSize',24);
set(h,'FontWeight','bold');
xlabel('Time (ms)');

% reduce size of data vector by down sampling to make CWT more efficient
t2 = [];
pData2 = [];
idx = 5; % sampling interval to keep
for i = 1:k
    if ~mod(i,idx)
        t2 = [t2 t(i)];
        pData2 = [pData2 pDataBarM(i)];
    end
    
end

%%% plot continuous wavelet transform %%%
subplot(3,1,[2 3]);
[wt,f] = cwt(pData2,Fs/idx);
[m,n] = size(wt);
h = pcolor((repmat(t2,m,1))*1000+tshift, repmat(f,1,n), abs(wt));
shading interp
colormap parula
ylim([5 220]);
set(h.Parent,'XTick',[0 100 200 300 400 500]);
xlim([0 500]);
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
h = text(0.006*1000, 208, '(b)');
set(h,'FontSize',24);
set(h,'Color',[1 1 1]);
set(h,'FontWeight','bold');

% plot dominant frequency
f0 = 11; % dominant frequency
for i = 1:3
    h = hline(i*f0);
    h.Color = [1 1 1];
    h.LineStyle = ':';
    h.LineWidth = 2;
end



