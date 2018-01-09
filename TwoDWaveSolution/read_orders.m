function [tel] = read_orders(telnum,order_data)
% Read in raw Zemax spectral format data according to Macro written on 9/29/16
% Parse and Reshape data to produce upright 2D arrays corresponding to detector
% Add conversion in the future

order_data_old = order_data;

%offsets measured by zooming into footprint diagram 12.2 11/1/17
%x offset = -10.578, y offset = -1.84

% order_data(:,2) = order_data(:,2)-10.578;
% order_data(:,3) = order_data(:,3)-1.84;  
% order_data(1:4:end,:)=[];

tel_w = reshape(order_data(telnum:3:end,1),[9,36])';
tel_y = reshape(-order_data(telnum:3:end,2),[9,36])';
tel_x = reshape(order_data(telnum:3:end,3),[9,36])';

% tel = zeros(39,9,3);

tel(:,:,1) = tel_x;
tel(:,:,2) = tel_y;
tel(:,:,3) = tel_w;

% Test figure
% figure(1)
% hold on
% plot(tel_x,tel_y,'ok')

end
