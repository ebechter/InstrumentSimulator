function [p1,ret,chebs] = WaveSolution(file,trace) 


fprintf('Fitting Zemax data with 2D Chebyshev\n') 
order_data=dlmread(file);

order_data=flipud(order_data);
polyorder = 4;


% choose which order to analyze  1,2,3

tel = fliplr(read_orders(trace,order_data));


% ----- Fit each order with a 2nd Order polynomial -----% 

% figure
% hold on
for ii = 1:36

 a{ii} = polyfit(tel(ii,:,1),tel(ii,:,3),polyorder);
%  plot(tel(ii,:,1),tel(ii,:,3),'.k')
%  plot(-25:25,polyval(a{ii},-25:25))
 
end
% box on
% grid on
% ylabel('Wavelength (\mum)')
% xlabel('X (pixels)')
% title('4th order polynomial fit')


for ii = 1:36
lam_cent(ii) = polyval(a{ii},0);
end
[order0] = find_m(lam_cent');

% Chebyshev fit parameters 
Inverse = 1;
npix = 4096;
ntotal = 39;
nx = 7;
nm = 7;
wavelengths = [];
pix_centers =[];
orders = [];
nparams = (min(nx,nm) + 1) .* (2*max(nx,nm) - min(nx,nm) + 2)./ 2;
p0 = zeros(1,nparams);

for ii = 1:36
wavelengths = [wavelengths tel(ii,:,3)];
pix_centers = [pix_centers tel(ii,:,1)];
orders = [orders (39+1-ii)*ones(size(tel(ii,:,1)))];
end


% ----- Fit zemax with 2D Chebyshev -----% 

[p1,ret, chebs] = Global_Wavelength_Solution(pix_centers, wavelengths, orders, p0,...
    order0, ntotal, npix, Inverse,nx,nm);
end




