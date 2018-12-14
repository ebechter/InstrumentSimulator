close all;
clear
addpath('/Users/ebechter/Desktop')

saving = 1;
% load data

% % Figure 1
% makePlot('Al',saving)
% 
% % Figure 2 
% makePlot('detH4RG',saving)

% Figure 3 
makePlot('Fathom',saving)

% Figure 4 
% makePlot('Xdisp',saving)




function [] = makePlot(prefix,saveflag)
S = dir([prefix,'*.mat']);
files = extractfield(S,'name');
color = viridis(length(files)+1);
figure
set(gcf, 'Position',  [1489          25        1225         400])
for i = 1:length(files)
    a = load(files{i},'-mat');
    fields = fieldnames(a);
    h{i}=plotOneCurve(a.(fields{1})(:,1),a.(fields{1})(:,2)/100,color(i,:),2,'Wavelength (nm)','Efficiency');
    legendnames{i} = files{i}(1:end-4);
end
legend(legendnames,'location','best')
ax = gca;
ax.LineWidth = 1.5;    
ax.FontSize = 20;
ylim([0 1])
xlim([a.(fields{1})(1,1) a.(fields{1})(end,1)])
grid minor
if saveflag
    tightfig(gcf)
    savefig(gcf,['figures/' prefix])
%     fig=gcf;
%     fig.PaperPositionMode = 'auto';
%     fig_pos = fig.PaperPosition;
%     fig.PaperSize = [fig_pos(3) fig_pos(4)];
%     print(fig, ['figures/' prefix],'-dpdf')
end
end


function [h] = plotOneCurve(x,y,color,linewidth,xlab,ylab)
hold on
h = plot(x,y,'-','color',color,'linewidth',linewidth);
box on
grid on
xlabel(xlab)
ylabel(ylab)
end