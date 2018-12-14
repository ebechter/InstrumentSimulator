% clear 
% close 

%Grating
path = 'S:\Simulator\RefFiles\Curves\RawCurves\echelle80k7_Littrow\';
file_id = strcat(path,'*.dat');
allFiles = dir(file_id);
fname = {allFiles.name};



% for ii = 1:1
% [orderData] = importGratingfile(fname{ii});
% [waveData] = orderData(:,1);
% [eff] = (orderData(:,2)+orderData(:,3))./2;
% end


waveData =[];
eff =[];

for ii = 1:size(fname,2)
[orderData] = importGratingfile(fname{ii});
if orderData(149,1) == orderData(150,1)
    orderData(149,1) = mean(orderData(148:150,1));
end
waveData = [waveData,orderData(:,1)];
eff = [eff,(orderData(:,2)+orderData(:,3))./2];
end

ub = max(max(waveData(waveData~=0)));
lb = min(min(waveData(waveData~=0)));
grating(:,1)= lb:1e-10:ub;

figure
hold on
for ii = 1:size(fname,2)
    temp = waveData(:,ii);
    ind = (temp ~=0);
    plot(waveData(ind,ii),eff(ind,ii))
    
    grating(:,ii+1) = interp1(waveData(ind,ii),eff(ind,ii),grating(:,1),'pchip',0);
end

return

%Grabit curves
filename = 'ICOS_ADC_AR.mat';
file = filename;
load(file)
curve = ICOS_ADC_AR;
x(:,1) = curve(:,1);
yr(:,1) = curve(:,2);
xq = 860:1:1545;
vq=interp1(x,yr,xq,'pchip','extrap');
x = permute(xq,[2,1]);
yr = permute(vq,[2,1]);
yt = (1-yr);

figure
% plot(x,yr)
% hold on
plot(x,yt)

mean(yt)

Out(:,1) = x;
Out(:,2) = yt;

ICOS_ADC = Out;
save('ICOS_ADC.mat','ICOS_ADC')

return

%Thorlabs curves
[output] = importThorlabsData('FB1330-12.xlsx');
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