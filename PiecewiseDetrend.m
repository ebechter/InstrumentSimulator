%% Function -- PiecewiseDetrend  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Normalize (flatten) the spectrum in preparation for cross correlation. This function removes large scale
% effects across an entire order after the removal of the blaze function. Samples of 
% the continuum are selected using a combination of a floor value to avoid some deep absorption lines or zero values and 
% the laplician of the flux to identify values that are not part of a steep slope in the spectrum (halfway down an 
% absorption line, for example). The samples are fit with a polynomial which is then divided out.  
%
% NOTES:
% 
% INPUT:
% spectrum - 2D array containing wavelength and flux values for a spectrum segment (e.g. one order)
% floorValue - scalar value used ot limit flux samples
% Del2Value - scalar value used to limit laplacian samples
% 
% OUTPUT:
% NormSpectrum - Continuum normalized spectrum 
%
% WRITTEN: Eric Bechter, 2016a.
% REVIEWED: 
% TESTED:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPROVEMENTS: (1) A very steep gradient across the entire order would require refining the floor value.
%               (2) Variable polynomial and scalar values? or as inputs to the function?
%               (3) If you don't get enough samples, loosen restrictions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectrum] = PiecewiseDetrend(spectrum)
% x = 3;
% xsamples = spectrum([1:x length(spectrum(:,1))-x:length(spectrum(:,1))],1);
% ysamples = spectrum([1:x length(spectrum(:,1))-x:length(spectrum(:,2))],2);
global colors


N=20;
split = [1, 20:floor(length(spectrum(:,2))/N):length(spectrum(:,2))-40, length(spectrum(:,2))-40,length(spectrum(:,2))];

var = [];
jj = 1;
for ii = 1:length(split)-1
    
    
    [~,ind(ii)] = max(spectrum(split(ii):split(ii+1),2));
    if ii>1
        ind(ii)=ind(ii)+split(ii)-1;
    end
end

ind = unique(ind);

for ii = 2:length(ind)-1
    if ind(ii)-ind(ii-1) <=100
        a=spectrum(ind(ii),2);
        b=spectrum(ind(ii-1),2);
        [~,test]=max([a,b]);
        
        if test == 1
             
            var(jj) = ii-1;
        else
            var(jj) = ii;
            
        end
        
        jj = jj+1;
    end
end
var=unique(var);
ind(var)=[];
% ind = [1,ind,length(ind)];
% ind = unique(ind);
xsamples = spectrum(ind,1);
ysamples = spectrum(ind,2);


% p = polyfit(xsamples,ysamples,1);  % This normalizes and starts to detrend
% trend = polyval(p,spectrum(:,1));

% trend = spline(xsamples,ysamples,spectrum(:,1));
trend = interp1(xsamples,ysamples,spectrum(:,1),'linear','extrap');


% figure()
% plot(spectrum(:,1),trend,'--k')
% hold on
% plot(spectrum(:,1),spectrum(:,2))
% plot(xsamples,ysamples,'o')
% hline = refline(0,1);
% hline.Color = 'k';
% 
spectrum(:,2) = spectrum(:,2)./trend;
% 
% figure
% plot(spectrum(:,1),spectrum(:,2),'color',colors{4})
% hline = refline(0,1);
% hline.Color = 'k';
% hline.LineStyle = '--';
end


% ind = find(flux>0.95);                    % doing this with only coninuum points helps to flatten
% 
% p2 = polyfit(wavelength(ind),flux(ind),1);
% line = polyval(p2,wavelength);
% flux = flux./line;

% 
% delta = max(flux)-1;
% 
% if max(flux) > 1
%     flux = flux-delta;
% end