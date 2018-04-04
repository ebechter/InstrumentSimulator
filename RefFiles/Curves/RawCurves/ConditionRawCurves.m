% clear 
% close 

%Thorlabs curves
[output]= importThorlabsData('FB1330-12.xlsx');
e = size(output,1);
FB1330_12(:,1) = flipud(output(2:e-2,3));
FB1330_12(:,2) = flipud(output(2:e-2,4)./100);

save('FB1330_12.mat','FB1330_12')

return

% clear 
% close 
% 
% %Entrance window from ECI
% filename = 'EntWindowTrans.mat';
% file = filename;
% load(file)
% curve = EntWindowTrans;
% x(:,1) = curve(:,1);
% yt(:,1) = curve(:,2)/100;
% xq = 401:1:1999;
% vq=interp1(x,yt,xq);
% x = permute(xq,[2,1]);
% yt = permute(vq,[2,1]);
% yr = (1-yt);
% 
% figure
% plot(x,yr)
% hold on
% plot(x,yt)
% 
% Out(:,1) = x;
% Out(:,2) = yt;
% 
% EntWTrans = Out;
% save('EntWTrans.mat','EntWTrans')

% clear 
% close 
% 
% %Fiber Dichroic
% filename = 'FiberDichroicRefl.mat';
% file = filename;
% load(file)
% curve = FiberDichroicRefl;
% x(:,1) = curve(:,1);
% yr(:,1) = curve(:,2)/100;
% xq = 401:1:1999;
% vq=interp1(x,yr,xq,'linear','extrap');
% x = permute(xq,[2,1]);
% yr = permute(vq,[2,1]);
% yt = (1-yr);
% 
% figure
% plot(x,yt)
% hold on
% plot(x,yr)
% 
% Out(:,1) = x;
% Out(:,2) = yr;
% 
% FiberDRefl = Out;
% save('FiberDRefl.mat','FiberDRefl')

clear
close 

%fathom gold file format:
% WL(nm),10deg r-pol,10deg s-pol,10deg p-pol,25deg r-pol,25deg s-pol,25deg p-pol
filename = 'S:\Simulator\RefFiles\Curves\FathomGold.mat';
file = strcat(filename);
load(file)

FathomGold10 = FathomGold(:,1:4);
FathomGold25 = FathomGold(:,[1,5:7]);

return
clear 
close 

%Fiber Dichroic
filename = 'QuadDichroicTrans.mat';
file = filename;
load(file)
curve = QuadDichroicTrans;
x(:,1) = curve(:,1);
yr(:,1) = curve(:,2)/100;
xq = 401:1:1999;
vq=interp1(x,yr,xq);
x = permute(xq,[2,1]);
yt = permute(vq,[2,1]);
yr = (1-yt);

figure
plot(x,yr)
hold on
plot(x,yt)

Out(:,1) = x;
Out(:,2) = yt;

QuadDTrans = Out;
save('QuadDTrans.mat','QuadDTrans')