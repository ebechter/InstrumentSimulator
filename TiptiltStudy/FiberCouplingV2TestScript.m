clear all
close all

% Fiber coupling class test script for V2 update. 

%========== Fiber Inputs ===========%

wfe(1) = 0;     % piston
wfe(2) = 0.4*sqrt(2);     % distortion/tilt
wfe(3) = 0;     % -
wfe(4) = 0.4*sqrt(3);  % defocus
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

% set fiber coupling parameters (used to vary coupling efficiency)

wavelength =flipud([1.29882 1.28752 1.27643 1.26552 1.25479 1.24425 1.23388 1.22368 1.21365 1.20378 1.19408 ...
    1.18452 1.17512 1.16587 1.15676 1.14779 1.13896 1.13027 1.12171 1.11327 1.10497 1.09678 1.08872 1.08077 ...
    1.07294 1.06522 1.05761 1.05011 1.04271 1.03542 1.02823 1.02114 1.01415 1.00725 1.00044 9.9373E-001 ...
    9.8710E-001 9.8057E-1 9.7411E-1]');
            
dispersion = false;

adc = 0;%60; % zenith angle for adc 0-60 in steps of 5

fiberpos = 5.8/40.6.*[0,0,0];% 1*[4/7,4/7,0]; % global position offset in microns (x,y,z)

dof = 0; % depth of focus (not sure if used yet)

[A] = FiberCouplingV2(wavelength,wfe,adc,fiberpos,dof,dispersion);

dz = wfe(4).*8.*wavelength.*(4^2);

figure
plot(wavelength,A.Rho./0.79)
