clear; clc;

% This is test number 2
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

%---------------------%
% Processing inputs
%---------------------%

parflag = true; %true or false

numworkers =2;

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

%sources{1}.throughput = {'lbt','lbti','fiberCh','SMFCoupling','fiberLink','spectrograph'};

sources{1}.throughput = {'lbt','demonstrator'};%'SMFCoupling','fiberLink','spectrograph'};

sources{1}.AO = 1;  

sources{1}.tput_flag = 1; % turn throughput on with 1

%---------------------%
% Observational inputs
%---------------------%

spType = 'M0V';

vmag = 0;

%V-I --> [F5 = 0.506, G0 = 0.664, G5 = 0.738,K0 = 0.853,K5 = 1.246, M0 = 1.848, M4=2.831, M7=4.52];

epsilon = 0;

vsini = 2.5; %km/s

rv = 0; %linspace(-30e3,30e3,20);

units = 'counts';

params.starInfo = {spType,vmag,epsilon,vsini,rv,units};

aoType = 'FLAO'; % 'FLAO' or 'SOUL'

entWindow = ''; %'ilocater' (Alluxa) or  empty for default lbti ZnSe []

zenith = 0*(pi/180); %radians

seeing = 1.1; %arcsec

conditions = {zenith,seeing,entWindow,aoType};

%--------------------%
% Persistence inputs
%--------------------%

persistentSource.name = 'star'; 

persistentSource.atmosphere = 0;

persistentSource.throughput = {'lbt','lbti','fiberCh','eta','fiberLink','spectrograph'};

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

SpecOrImager = 'Imager' ; %'Imager' or Spectrograph

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

optical.bandPass = 965:1285; % typical bandpass for fiber coupling

optical.wfef = [0,0,0,0,0,0,0,0]; % static aberrations at the fiber tip (same format as above)

optical.adc = 0; % zenith angle for adc 0-60 in steps of 5

optical.fiberpos = 1*[4/7,4/7,0]; % global position offset in microns (x,y,z)

optical.dof = 0; % depth of focus (not sure if used yet)

%---------------------%
% Output options
%---------------------%

pathprefix = pwd;

mag = 1.272;
for ii = 1:length(mag)
    
    vmag = mag(ii)+1.848; %M4 type
    
    starInfo = {spType,vmag,epsilon,vsini,rv,units};
    
    fname = [spType,'_',num2str(vmag)];
    
    fitsname = [pathprefix,'/Output/FiberCouplingV2',fname,'.fits'];
    
    headerinfo = {
    'version',version,'version';...
    'footprint',footprint,'spectrograph fpt';...
    'parallel',parflag,' ';...
    'scale',scale,'upsample factor';...
    'nOrders',nOrders,'num orders';...
    'aoType',aoType,'FLAO or SOUL AO';...
    'tracenum',tracenum','traces';...
    'entWindow',entWindow,'entrance window';...
    'zenang',zenith,' zenith in rad';...
    'zenseeing',seeing,'seeing in arcsec';...
    'adcpos',optical.adc,'ADC pos in deg'
    'fiberpos',optical.fiberpos,'position in mas';...
    'Inst',SpecOrImager,'instrument type';...
    'SpType',spType,'spec type';...
    'units',units,'units of spectrum';...
    'Vmag', vmag,' ';...
    'RV', rv,'Injected RV (m/s)';...
    'Vsini', vsini,' ';...
    'Tlrcs', sources{1}.atmosphere ,'Are tellurics included 1/0';...
    };
    
    [T] = SimulationMainThroughput(parallelInfo,runInfo,specInfo,fitsname, ...
    optical,starInfo,sources,SpecOrImager,conditions,headerinfo,persistentSource);
    
end

