%Ruilier Plot
function [] =  RuilierPlot(fc,a,colors)

figure()
box on
grid on
ylim ([0,1])
xlim([0,3])
ylabel('Normalized Coupling Efficiency')
xlabel('\sigma_{rms} (rads)')
hold on
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 16;

b=[1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8),1/sqrt(8),1/sqrt(5),1/sqrt(8),1/sqrt(8),1/sqrt(10),1/sqrt(10),1/sqrt(12),1/sqrt(12),1/sqrt(7)];

for kk = 2:10
    h(kk) = plot(a*b(kk),fc(kk,:),'linewidth', 1.5);
end

%--------%
% Strehl
%--------%

a2 = linspace(min(a),max(a),100);
h(1) = plot(a2,exp(-a2.^2),'linewidth', 1);
h(1).LineStyle = 'none';
h(1).Marker = 'd';
h(1).Color = colors{1};
h(1).MarkerFaceColor = colors{1};
h(1).MarkerSize = 4;

%--------%
% Tip-tilt
%--------%

h(2).LineStyle = ':';
h(2).Color = colors{2};
h(2).LineWidth = 1.5;
h(3).LineStyle = 'none';

%--------%
% Focus
%--------%

h(4).LineStyle = '--';
h(4).Color = colors{3};

%-------------%
% Astigmatism
%-------------%

h(5).LineStyle = '-.';
% h(5).Marker = 'x';
h(5).MarkerSize = 8;
h(5).Color = colors{4};
h(6).LineStyle = 'none';


%-------------%
% Coma
%-------------%

h(7).LineStyle = '-.';
h(7).Color = colors{5};
h(8).LineStyle = 'none';

%-------------------------%
% Spheric
%-------------------------%

% h(9).Marker = 'o';
h(9).Color = colors{6};

%-------------------------%
% Trefoil (triangular coma)
%-------------------------%

% h(10).Marker = '*';
h(10).Color = colors{7};
% h(11).LineStyle = 'none';

%-------------------------%
% Others
%-------------------------%
% h(12).Marker = 'p';
% h(13).Marker = 'p';
% h(14).LineStyle = '-.';
% h(15).LineStyle = '-.';
% h(16).Marker = 's';

legend([h(1),h(2),h(4),h(5),h(7),h(10),h(9)],...
    'strehl','tip/tilt (Z2)','defocus (Z4)','astigmatism (Z5)','coma (Z7)','trefoil (Z9)','spheric (Z11)')

  %'trefoil','trefoil','secondary astig','secondary astig','secondary coma','secondary coma','secondary spherical')
  
end
