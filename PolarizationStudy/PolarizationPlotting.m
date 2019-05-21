%Polarization Grating Plots
clear
%%------------
% Load Files
%%------------

load('S:\Simulator\PolarizationStudy\PolarizationTests\DichroismSet1.mat')

% -------------%
% Index Shifting
% -------------%
n = length(spectral_cell); % files
jj = 15; %order number;
kk = 1:n;
ind = kk(1:2:end);% index of every nth plot

% -----------%
% Plotting 
% -----------%

quart = char(188);     % 1/4 symbol
half = char(189);      % 1/2 symbol
quart3 = char(190);    % 3/4 aymbol
piStr = char(960);     % Greek letter pi

N = length(ind);

c = viridis(N);
c = c(1:N-2,:);
% c = flipud(c);

% -------------------------------%
% Grating Wavelngth vs. Intensity
% -------------------------------%

figure()
hold on

for ii= 1:N
    kk = ind(ii);
    h(ii) = plot(OrderWave{kk}{jj}(1,:),OrderFlux{kk}{jj}(1,:));
    D = Spectrograph{kk}.polarization(1,2);
    labels{ii} = ['$S = ',num2str(D,'%.2f'),'$'];
end



%settings
opt = [];
opt.MarkerSpacing = [100,100];
opt.LineWidth = 1.2*ones(N,1);
opt.Markers = {'none','none','none','none'};
opt.MarkerSize = [5,7];
opt.LineStyle = {'-','-','-','-'};
opt.Colors = c; % change plot color
opt.XLabel = '$Wavelength, \lambda, (nm)$'; % xlabel
opt.YLabel = '$Intensity , (cts)$'; %ylabel
% opt.XLim = [0.97,0.978];

% create the plot
setPlotProp(opt);
ax = gca;

[l,icons,~,~] = legend(labels,'Interpreter','latex','location','nw');
l.Box = 'off';
l.FontSize = 12;

% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'Figures/Dichroism','-dpdf')





