


file = '/Volumes/Software/SimulatorV5/RefFiles/Zemax/9.88_Format_jc.txt';



%% Read in spectral format file.
order_data=dlmread(file);
order_data=flipud(order_data);
% 1 Telescope data at a time
tel1 = read_orders(1,order_data); 
tel2 = read_orders(2,order_data);
tel3 = read_orders(3,order_data);

% telA follows this: (:,:,ii), ii = x,y,lambda

% Nanometers and bottom left 0 instead of microns and 

%% XY
% figure()
% hold on





for ii = 1:39
    %     hold on
    %
    %     % fit each order with a polynomial, save the coefficients
    yx_coeff(ii,1:3,1) = polyfit(tel1(ii,:,1),tel1(ii,:,2),2);
    yx_coeff(ii,1:3,2) = polyfit(tel2(ii,:,1),tel2(ii,:,2),2);
    yx_coeff(ii,1:3,3) = polyfit(tel3(ii,:,1),tel3(ii,:,2),2);
    
    % %       Plot all the orders
    %     plot(tel1(ii,:,1),polyval(xy_coeff(ii,1:3,1),tel1(ii,:,1)))
    %     plot(tel2(ii,:,1),polyval(xy_coeff(ii,1:3,2),tel2(ii,:,1)))
    %     plot(tel3(ii,:,1),polyval(xy_coeff(ii,1:3,3),tel3(ii,:,1)))
    %     pause
end
%% lam-X
% figure()
for ii = 1:39
    %     hold on
    %     plot(tel1(ii,:,1),tel1(ii,:,3),'o')
    
    % fit each order with a polynomial, save the coefficients
    lamx_coeff(ii,:,1) = polyfit(tel1(ii,:,1),tel1(ii,:,3),2);
    lamx_coeff(ii,:,2) = polyfit(tel2(ii,:,1),tel2(ii,:,3),2);
    lamx_coeff(ii,:,3) = polyfit(tel3(ii,:,1),tel3(ii,:,3),2);
    %
    %       Plot all the orders
    %     plot(tel1(ii,:,1),polyval(lamx_coeff(ii,:,1),tel1(ii,:,1)))
    %     plot(tel2(ii,:,1),polyval(lamx_coeff(ii,:,2),tel2(ii,:,1)))
    %     plot(tel3(ii,:,1),polyval(lamx_coeff(ii,:,3),tel3(ii,:,1)))
end
%% X-lam -- using 3rd order polynomial
% figure()
% for ii = 1:39
%     %     hold on
%     %     plot(tel1(ii,:,3),tel1(ii,:,1),'o')
%     
%     % fit each order with a polynomial, save the coefficients
%     xlam_coeff(ii,:,1) = polyfit(tel1(ii,:,3),tel1(ii,:,1),2);
%     xlam_coeff(ii,:,2) = polyfit(tel2(ii,:,3),tel2(ii,:,1),2);
%     xlam_coeff(ii,:,3) = polyfit(tel3(ii,:,3),tel3(ii,:,1),2);
%     %
%     % %       Plot all the orders
%     %     plot(tel1(ii,:,3),polyval(xlam_coeff(ii,:,1),tel1(ii,:,3)))
%     %     plot(tel2(ii,:,3),polyval(xlam_coeff(ii,1:3,2),tel2(ii,:,3)))
%     %     plot(tel3(ii,:,3),polyval(xlam_coeff(ii,1:3,3),tel3(ii,:,3)))
% end



tel_data(1,:,:,:) = tel1;
tel_data(2,:,:,:) = tel2;
tel_data(3,:,:,:) = tel3;




lx = lamx_coeff(1,:,3);

x = (-20:20);
lam = polyval(lx,x);

figure()
plot(diff(lam),'.');

% 
% lam = (970:980)*1e-3;

a = lx(1);
b = lx(2);
c = lx(3);

xl = (-b + sqrt(b^2 - 4*a*(c-lam)))/(2*a);







