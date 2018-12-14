a = dir('*.fig');
fignames = extractfield(a,'name');
for i = 1:length(fignames)
    openfig(fignames{i})
    fig=gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,[fignames{i}(1:end-4) '.pdf'],'-dpdf')
    close all
    
end