%Ruilier Study 
%Purpose; Use the fiber coupling class to reproduce and expand on the work Ruilier et al
%2001

clear; clc; close all;

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

%---------------------%
% Wavefront error
%---------------------%

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

%---------------------%
% Fiber Inputs
%---------------------%

% set fiber coupling parameters (used to vary coupling efficiency)

wavelength =[1.02823]'; %microns

dispersion = false;

adc = 0; % zenith angle for adc 0-60 in steps of 5

fiberpos = 0*[1.4,1.4,0];% 1*[4/7,4/7,0]; % global position offset in microns (x,y,z)

dof = 0; % depth of focus (not sure if used yet)

a = [0:0.1:10]; %0.5~0.7 microns %a=1 in tip and tilt account for 30%

for jj = 2%length(wfe)
    
    f =[];
    rmsWFE = [];
    rmsPhase = [];
    X = [];
    Y= [];
    for ii = 1:length(a)
        
        wfe(jj) = a(ii);
        
        [A] = FiberCouplingV2(wavelength,wfe,adc,fiberpos,dof,dispersion);
        
        [x,y] = centroid(A.PSF.*conj(A.PSF));
        
        xpos(ii) = x*mean(mean(diff(A.FPgridx')));
        
        ypos(ii) = y*mean(mean(diff(A.FPgridy)));
        
        f = [f; A.Rho];
        
        rmsPhase= [rmsPhase; A.rmsPhase];
        
        rmsWFE = [rmsWFE; A.rmsWFE];
        
        phi(ii)=2*[(a(ii)./1.064e-6)-1i]*pi;
        
        X = [X,x];
        Y = [Y,y];
        
    end
    
    return
    a2 =[0:0.01:max(a)];
    
    f2 = interp1(a,f,a2);
    
    wfe(jj)=0;
    
    fc(jj,:) = f2; % for multiple wavelengths
end

fc= fc./max(fc);

save('fiber coupling', 'fc','a2')


%---------------------%
% Output Plots
%---------------------%

science = viridis(8);
for ii = 2:8
    colors{ii} = science(ii,:);
end

RuilierPlot(fc,a2,colors)

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


