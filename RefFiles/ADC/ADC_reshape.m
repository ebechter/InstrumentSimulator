% ADC residual dispersion and shifting centroids
clear all
close all

%% Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end
colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};
colors{11}= colors{3};
colors{12}= colors{3};
colors{13}= colors{3};
colors{14}= colors{3};

clear d

% rows sort by wavelength wavlengths
% layers are different airmasses
% columns are fields (wav,x,y)

load('S:\SimulatorV10\RefFiles\ADC\ADC.mat')
[r,c] = size(ADC);
nlay = 13;
adc = permute(reshape(ADC',[c,r/nlay,nlay]),[2,1,3]);
% po = zeros(12,8);
pscale = 1;%6.9;% mas/mum
for ii = 1:13
    
%     adcn = adc(:,3,:)-mean(adc(:,3,:));
%        
%     figure(1)
%     hold on
%     
%     if ii == 12
%         p0{ii} = polyfit(adc(:,1,ii),adcn(:,1,ii)*1e3*pscale,7);
%     else
%         p0{ii} = polyfit(adc(:,1,ii),adcn(:,1,ii)*1e3*pscale,2);
%         
%     end
    
    %mean centered
    adcn(:,1,:) = adc(:,3,:)-mean(adc(:,3,:));
    adcn(:,2,:) = adc(:,2,:)-mean(adc(:,2,:));
    r = sqrt((adcn(:,2,:).^2 + adcn(:,1,:).^2));
    
    %short wavelength centered
%     adcn(:,1,:) = adc(:,3,:)-(adc(1,3,:));
%     adcn(:,2,:) = adc(:,2,:)-(adc(1,2,:));
%     r = sqrt((adcn(:,2,:).^2 + adcn(:,1,:).^2));
    
    figure(2)
    hold on
%     for ii = 1:13
    plot(adc(:,1,ii),r(:,1,ii)*1e3*pscale,'.-','color','k')
%     end

    p0{ii} = polyfit(adc(:,1,ii),r(:,1,ii)*1e3*pscale,2);
    x1 = linspace(0.970, 1.310, 20);
    y1 = polyval(p0{ii},x1);
    
%     plot(adc(:,1,ii),adcn(:,1,ii)*1e3*pscale,'o','color','k')
    if ii > 6
    a(ii) = plot(x1,y1,'linestyle','-.','color',colors{ii-6});
    else
    a(ii) = plot(x1,y1,'color',colors{ii});
    end
    hold off
end
box on
% hline = refline(0,4);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
% hline = refline(0,-4);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;

l = legend(a,{'0 deg','5 deg','10 deg','15 deg','20 deg','25 deg','30 deg','35 deg','40 deg','45 deg','50 deg','55 deg','60 deg'});
xlabel ('Wavelength (\mum)')
ylabel ('\Delta x (mas)')
l.FontSize = 14;

set(gca,'FontSize',16)

% ylim([-5,5])



% output 



