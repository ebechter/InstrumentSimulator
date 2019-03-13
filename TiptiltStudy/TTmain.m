%Simulated Fiber coupling time series

%Call the fiber coupling class

%Set the wavelength to monochromatic

%Inject the correct tip/tilt WFE
%Units should be in mas I think then convert to waves for phase.

% Test:

%Repeat according to a distribution or time series

%Output a fiber coupling efficiency signal

%Normalize and compare to ruilier tip/tilt losses

clear; clc; %close all

% Start from the correct path
current_path = pwd;
if strcmp(current_path(end-8:end),'Simulator')~=1
    
    if strcmp(current_path(2:8),'Volumes')==1
        cd '/Volumes/Software/Simulator/'
        
    elseif strcmp(current_path(2:4),'afs')==1
        cd '/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/'
        
    else
        cd([current_path(1:2) '\Simulator'])
    end
end
clear current_path

addpath(genpath(pwd))
tic

%========== Wavefront error ===========%

% Wavefront error generation using Zernike % indexing starts at 0 for Zernike #s. 0 is piston... etc. follows Wyant scheme

wfe = zeros(16,1); % static aberrations in the spectrograph

% characteristic wavefront error in radians describing the amplitude (not RMS or P-V)

wfe(1) = 0;     % piston
wfe(2) = 0;     % distortion/tilt
wfe(3) = 0;     % -
wfe(4) = 0;     % defocus
wfe(5) = 0;     % Primary astigmatism
wfe(6) = 0;     % -
wfe(7) = 0;     % Primary coma
wfe(8) = 0;     % -
wfe(9) = 0;     % Spherical abb
wfe(10) = 0;    % Elliptical coma
wfe(11) = 0;    % -
wfe(12) = 0;    % Secondary astigmatism
wfe(13) = 0;    % -
wfe(14) = 0;    % Secondary coma
wfe(15) = 0;    % -
wfe(16) = 0;    % Secondary spherical abb


%========== Fiber Inputs ===========%

% set fiber coupling parameters (used to vary coupling efficiency)

wavelength =[1.02823]'; %microns

dispersion = false;

adc = 0; % zenith angle for adc 0-60 in steps of 5

fiberpos = [0,0,0];% 1*[4/7,4/7,0]; % global position offset in microns (x,y,z)

dof = 0; % depth of focus (not sure if used yet)

a = 0* 5.8/40.6.*(-30:1:30); %0.5~0.7 microns %a=1 in tip and tilt account for 30%

for jj = 1:1
    
    f =[];
    
    for ii = 1:length(a)
        
        fiberpos(1) = a(ii);
        
        [A] = FiberCouplingV2(wavelength,wfe,adc,fiberpos,dof,dispersion);
        
        [x,y] = centroid(A.PSF.*conj(A.PSF));
        
        xpos(ii) = x*mean(mean(diff(A.FPgridx')));
        
        ypos(ii) = y*mean(mean(diff(A.FPgridy)));
        
        f = [f; A.Rho./0.79];
        
    end
    
    wfe(jj)=0;
    
    fc(jj,:) = mean(f,2);
end

xpos = xpos-xpos(1); %0 centered centroid
ypos = ypos-ypos(1); %0 centered centroid
% xpos = xpos*40.6/5.8e-6; %mas
% ypos = ypos*40.6/5.8e-6; %mas

amp = 40.6/5.8.*a;
amph =(min(amp):0.01:max(amp));

% amp = a;
% amph =[0:0.01:max(amp)];
fch = interp1(amp,fc,amph);

figure
plot(amph,fch)
hold on
plot(amp,fc)
 %symmetric offset if centroid data correlates. Misalignment would be a mix
 %kind of. 
 
return
figure
hold on
b=[1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8)];
for kk = 1:2
    plot(a*b(kk),fc(kk,:))
end
legend('tip','tilt','defocus','astig','astig','coma')


function [xc, yc] = centroid(img)
% CENTROID computes image centroid location
%
% [xc, yc] = CENTROID(img)
%
% Parameters
% ----------
% img : m x n array
%   Image
%
% Returns
% -------
% x : float
%   x centroid (column)
%
% y : float
%   y centroid (row)
%

[rows, cols] = size(img);

xc = floor(cols/2) + 1;
yc = floor(rows/2) + 1;

m  = sum(sum(img));

if m ~= 0
    mx = sum(img)  * (1:cols)';
    my = sum(img') * (1:rows)';
    
    xc = mx / m;
    yc = my / m;
end
end


