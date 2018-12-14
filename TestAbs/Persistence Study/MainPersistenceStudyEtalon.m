% simulation Main Script
% clear up the workspace
clear; clc; close all
addpath(genpath([pwd '/../']))
tic

% Parameters 
version = '11';
footprint = '12.22';
numworkers =4;
parflag = false; %true or false
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
rv = 0;


units = 'counts';


persistentRV = 0;%[0:500:5000 7000:2000:1e4 2e4:1e4:3e4];
persistentRV = [-fliplr(persistentRV(2:end)) persistentRV]; 


persistentSource.name = 'etalon'; 
persistentSource.atmosphere = 0;
persistentSource.throughput = {};
persistentSource.AO = 0;  
persistentSource.spType = '';
persistentSource.vmag = 11.8;
persistentSource.epsilon = 0;
persistentSource.vsini = 1;
persistentSource.units = 'counts';
persistentSource.persistence = 0.1; 

% 
% % Source 2  
% % ----------- % 
% sources{2}.name = 'etalon'; 
% sources{2}.atmosphere = 0;
% sources{2}.throughput = {'calibration','spectrograph'};
% sources{2}.AO = 0; 


% Source 3  
% ----------- % 
% sources{3}.name = 'flat'; 
% sources{3}.atmosphere = 0;
% sources{3}.throughput = {'calibration','spectrograph'};
% sources{3}.AO = 0; 

pathprefix = pwd;

for ii = 1:length(persistentRV)

parallelInfo = {numworkers,parflag};
runInfo = {version,scale};
conditions = {zenith,seeing,entWindow,aoType};
specInfo = {footprint,nOrders,cheby,p1,ret,tracenum,order_coeff,wave_coeff};
starInfo = {spType,vmag,epsilon,vsini,rv,units};

persistentSource.rv = persistentRV(ii);

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
    'RV', rv,'Injected RV (m/s)';...
    'Vsini', vsini,' ';...
    'Tlrcs', sources{1}.atmosphere ,'Are tellurics included 1/0';...
    'Pers',persistentSource.persistence,'level of persistence';...
    'PRV', persistentSource.rv,'persistent source RV';...
    'PSpType',persistentSource.spType,'persistent source type';...
    'Pvmag',persistentSource.vmag,'persistent source mag' ;...
    };


fname = ['PersistentEtalon_0.1_' num2str(ii)];
fitsname = [pathprefix,'/../Output/Persistence/',fname,'.fits'];

%Call the main script

SimulationMainP(parallelInfo,runInfo,specInfo,fitsname, ...
    wfe,starInfo,sources,polarization,SpecOrImager,conditions, headerinfo,persistentSource)
toc
end


