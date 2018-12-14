clear; clc; %close all
addpath(genpath(pwd))
tic

%---------------------%
% Processing inputs
%---------------------%

parflag = true; %true or false

numworkers =4;

parallelInfo = {numworkers,parflag};

version = '11';

scale = 1;

footprint = '12.22';

runInfo = {version,scale};

%---------------------%
% Source inputs
%---------------------%

sources{1}.name = 'star'; 

sources{1}.atmosphere = 1;

sources{1}.throughput = {'lbt','lbti','fiberCh'};

sources{1}.AO = 1;  

sources{1}.tput_flag = 0; % turn throughput on with 1

%---------------------%
% Observational inputs
%---------------------%

spType = 'M0V';

vmag = 10.8;

epsilon = 0;

vsini = 2.5; %km/s

rv = 0; %linspace(-30e3,30e3,20);

units = 'counts';

params.starInfo = {spType,vmag,epsilon,vsini,rv,units};

aoType = 'SOUL'; % 'FLAO' or 'SOUL'

entWindow = 'ilocater'; %'ilocater' (ECI) or  empty for default lbti []

zenith = 30*(pi/180); %radians

seeing = 1.0; %arcsec

conditions = {zenith,seeing,entWindow,aoType};

%--------------------%
% Persistence inputs
%--------------------%

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

%------------------%
% Optical inputs
%------------------%

% sets the instrument type, static wavefront error, poalrization and fiber
% parameters

%========== Instrument type ===========%

% Choose between simulating and imager or a spectrograph

SpecOrImager = 'Spectrograph' ; %'Imager' or Spectrograph

%========== Spectrograph info ===========%

nOrders = 36;

tracenum = 1;

cheby=0;

load polycoeffs2

load chebycoeffs2

specInfo = {footprint,nOrders,cheby,p1,ret,tracenum,order_coeff,wave_coeff};

%========== Wavefront error ===========%

% Wavefront error generation using Zernike % indexing starts at 0 for Zernike #s. 0 is piston... etc. follows Wyant scheme 

optical.wfe = zeros(16,1); % static aberrations in the spectrograph

% characteristic wavefront error in radians describing the amplitude (not RMS or P-V)

optical.wfe(1) = 0;     % piston
optical.wfe(2) = 0;     % distortion/tilt
optical.wfe(3) = 0;     % -
optical.wfe(4) = 0;     % defocus
optical.wfe(5) = 0;     % Primary astigmatism
optical.wfe(6) = 0;     % -
optical.wfe(7) = 0;     % Primary coma
optical.wfe(8) = 0;     % -
optical.wfe(9) = 0;     % Spherical abb
optical.wfe(10) = 0;    % Elliptical coma
optical.wfe(11) = 0;    % -
optical.wfe(12) = 0;    % Secondary astigmatism
optical.wfe(13) = 0;    % -
optical.wfe(14) = 0;    % Secondary coma
optical.wfe(15) = 0;    % -
optical.wfe(16) = 0;    % Secondary spherical abb

%========== Polarization State ===========%

% Can be used on a single optic which is set in instrument
% params model! This polarization parameter sets the inital state going into the
% system but is not updated after hitting the first polarization sensitive
% optic. Once the state is set to update the polarization propagation
% should work through the whole system. The grating is still sythetic but
% based on measured alumnium curves. S and P polarizations 

optical.polarization = [1,0.5,0]; % [degree of polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)]

%========== Fiber Inputs ===========%

% set fiber coupling parameters (used to vary coupling efficiency)

optical.bandPass = 965:1300; % typical bandpass for fiber coupling

optical.wfef = [0,0,0,0,0,0,0,0]; % static aberrations at the fiber tip (same format as above)

optical.adc = 30; % zenith angle for adc 0-60 in steps of 5

optical.fiberpos = 0*[4/7,4/7,13]; % global position offset in microns (x,y,z)

optical.dof = 0; % depth of focus (not sure if used yet)

%---------------------%
% Output options
%---------------------%

pathprefix = pwd;

for ii = 1:1
    
    starInfo = {spType,vmag,epsilon,vsini,rv(ii),units};
    
    fname = ['M0'];
    
    fitsname = [pathprefix,'/Output/PhNoiseFDR_Tellurics_OFF',fname,'.fits'];
    
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
    
    SimulationMain2(parallelInfo,runInfo,specInfo,fitsname, ...
    optical,starInfo,sources,SpecOrImager,conditions,headerinfo,persistentSource)
    
    save params
end

