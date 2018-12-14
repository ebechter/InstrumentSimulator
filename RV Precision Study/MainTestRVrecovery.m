% simulation Main Script
% clear up the workspace
clear; clc; %close all
addpath(genpath(pwd))
tic

% Parameters 
version = '11';
footprint = '12.22';
numworkers =4;
parflag = true; %true or false
scale = 1;
load polycoeffs2
load chebycoeffs2
nOrders = 36;
cheby=0;
aoType = 'SOUL'; % 'FLAO' or 'SOUL'
entWindow = 'ilocater'; %'ilocater' (ECI) or  empty for default lbti []
zenith = 10*(pi/180); %rad
seeing = 1.1; %arcsec
SpecOrImager = 'Spectrograph' ; %'Imager' or Spectrograph
tracenum = 1;


% Polarization state
% right now this can be used on a single optic which is set in instrument
% optical model! This polarization parameter sets the inital state going into the
% system but is not updated after hitting the first polarization sensitive
% optic. Once the state is set to update the polarization propagation
% should work through the whole system. The grating is still sythetic but
% based on measured alumnium curves. S and P polarizations 

polarization = [1,0.5,0]; % [degree of polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)]

% Wavefront error generation using Zernike % indexing starts at 0 for Zernike #s. 0 is piston... etc. follows Wyant scheme 

%--------------%
wfe = zeros(16,1);


% characteristic wavefront error in radians describing the amplitude (not RMS or P-V)
wfe(1) = 0;     % piston
wfe(2) = 0;     % distortion/tilt
wfe(3) = 0;     % -
wfe(4) = 0;     % defocus
wfe(5) = 0;     % Primary astigmatism
wfe(6) = 0;     % -
wfe(7) = 0;    % Primary coma
wfe(8) = 0;    % -
wfe(9) = 0;     % Spherical abb
wfe(10) = 0;    % Elliptical coma
wfe(11) = 0;    % -
wfe(12) = 0;    % Secondary astigmatism
wfe(13) = 0;    % -
wfe(14) = 0;    % Secondary coma
wfe(15) = 0;    % -
wfe(16) = 0;    % Secondary spherical abb


% Source Generation and Options 
%--------------%
% Can take up to 3 sources if you are doing spectroscopy, only 1 source for imaging
%
% source list
% 'star', 'etalon','flat','cal'
%
% Atmosphere flag: 
% 1 or 0
%
% throughput instrument list:
% 'lbt','lbti','fiberCh','imageCh','quadCh','wfc','fiberLink','spectrograph','calibration','filter'
%
% AO flag
% 1 or 0 
%
% Example: 
% sourceX = 'star'; 
% atmosphereX = 0;
% throughputX = {'lbt','lbti','fiberCh','fiberLink','spectrograph'};
% useAO1 = 1;
%--------------%

% Source 1  
% ----------- % 
sources{1}.name = 'star'; 
sources{1}.atmosphere = 0;
sources{1}.throughput = {'lbt','spectrograph'};
sources{1}.AO = 0;  
sources{1}.tput_flag = 0;

spType = 'M0V';
vmag = 11.8;
epsilon = 0;
vsini = 1 ;
rv = 0;%linspace(-30e3,30e3,20);
units = 'counts';

persistentSource.name = 'star'; 
persistentSource.atmosphere = 0;
persistentSource.throughput = {'lbt','lbti','fiberCh','fiberLink','spectrograph'};
persistentSource.AO = 1;  
persistentSource.spType = 'G2V';
persistentSource.vmag = 11.8;
persistentSource.epsilon = 0;
persistentSource.vsini = 1;
persistentSource.units = 'counts';
persistentSource.persistence = 0; 
persistentSource.rv =0;

parallelInfo = {numworkers,parflag};
runInfo = {version,scale};
conditions = {zenith,seeing,entWindow,aoType};
specInfo = {footprint,nOrders,cheby,p1,ret,tracenum,order_coeff,wave_coeff};


pathprefix = pwd;

for ii = 1:length(rv)
    
    starInfo = {spType,vmag,epsilon,vsini,rv(ii),units};

headerinfo = {
    'version',version,'version';...
    'footprint',footprint,'spectrograph fpt';...
    'parallel',parflag,' ';...
    'scale',scale,'upsample factor';...
    'nOrders',nOrders,'num orders';...
    'aoType',aoType,'FLAO or SOUL AO';...
    'tracenum',tracenum','traces';...
    'entWindow',entWindow,'entrance window';...
    'zenith',zenith,' zenith in rad';...
    'seeing',seeing,'seeing in arcsec';...
    'Inst',SpecOrImager,'instrument type';...
    'SpType',spType,'spec type';...
    'units',units,'units of spectrum';...
    'Vmag', vmag,' ';...
    'RV', rv(ii),'Injected RV (m/s)';...
    'Vsini', vsini,' ';...
    'Tlrcs', sources{1}.atmosphere ,'Are tellurics included 1/0';...
    };


fname = ['M0' num2str(ii)];
fitsname = [pathprefix,'/Output/',fname,'.fits'];
%Call the main script

SimulationMain(parallelInfo,runInfo,specInfo,fitsname, ...
    wfe,starInfo,sources,polarization,SpecOrImager,conditions, headerinfo,persistentSource)

end






