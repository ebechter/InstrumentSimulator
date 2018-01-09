% SpectrographThroughput

%% Global variables
addpath(genpath(pwd))
global allardpath telluricfile skybackfile starfile calibfile
global curve_dir
global colors
global strehl_dir
global zmx_kernel
global zmx_detector
global gratingfile
global zernikefile

%% Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end
colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};
clear d
%% Set paths (works on windows and mac)
pathprefix = pwd;
allardpath = [pathprefix '/../Spectral_Catalogs/FAllard/CIFIST6b_trimmed/'];
telluricfile =[pathprefix '/RefFiles/Atmosphere/telluric_200.mat'];
skybackfile = [pathprefix '/RefFiles/SkyBackground/SkyBackground.mat'];
% starfile =  [pwd '/RefFiles/Star/fullpecautmamajek.xlsx'];
starfile =  [pathprefix '/RefFiles/Star/fullpecautmamajek.xlsx'];
curve_dir = [pathprefix '/RefFiles/Curves/'];
strehl_dir= [pathprefix '/RefFiles/AO/'];
zmx_kernel  = [pathprefix '/RefFiles/Zemax/Kernel/PSF_256_025.txt'];
zmx_detector = [pathprefix '/RefFiles/Zemax/9.88_Format_jc.txt'];
gratingfile = [pathprefix '/RefFiles/Curves/GratingEff_highsamp.mat'];
calibfile = [pathprefix '/Calibration/Lorentzian.mat'];
zernikefile = [pathprefix '/RefFiles/AO/zernike_index.mat'];
%end of path block


band = (950:1350);

%Initialise dedicated object for calibration fiber
Spec = Calibration(band);
Spec = CustomizeBandpass (Spec,[950 1350]);%Customize Bandpass (all have defaults preloaded)
Spec = Trim_Throughput(Spec);%Trim to custom Bandpass (1,:)
Spec = Integrated_Transmission(Spec);
Spec = Include_Grating(Spec);%Include grating into throughput object by creating a new property Calibration.Orders

