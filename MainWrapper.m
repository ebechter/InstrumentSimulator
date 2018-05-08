


% Create extended color list
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end

colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};
colors{11} =[0 0.3 0];
clear d

%-----------------------Parallel Settings----------------------%
numworkers=4;
parflag = true;
poolobj = gcp('nocreate');

if parflag ==true  % 
    if isempty(poolobj)
        %no pool exists and I want one
        parpool(numworkers);
    else
        %just run it b/c a pool already exists
    end
else
    delete(poolobj); % should shut down the current pool or does nothing if there is no pool
end


%-----------------------User Inputs-----------------------------%
wavesolution = 0; % 1 = recalculate, 0 = load previous
cheby = 0; % 1 = use chebyshev, 0 = 4th order polynomial
norders=1; % max number of orders now 36
tracenum=[1]; %  [1,2,3] traces: 1 is (bot) telescope fiber, 2 (top) telescope fiber, 3 is cal fiber (middle). 

%-----Instrument Inputs-----% 
StarlightPath = 'Spectrograph'; % Options are "Fiber,Image,Quad,Spectrograph,WFC"
Throughput =1; %0 = No, 1 = Yes %includes all throughout effects (except for dispersion)
mode =3; %0 = Science (S,E,S), 1 = Flat (F,F,F) , 2 = (E,E,E), 3 = Testing (S, E, F); 
band = (900:1350);
entwindow = 'ECI'; %options are: 'ECI','sim','existing'
polarization = [0.5,0.1,0]; % degree of polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)

%-----FiberCoupling Inputs-----% 
wfe = [0,0,0,0,0,0,0,0,0,0,0]; % Pupil Zernike amplitudes piston, xtilt, ytilt, etc...
adc = 45; %Choose values from 0-60 in steps of 5 for the zenith angle
%The ADC will automatically incorporates the residual disperion into the beam based on the airmass
fiberpos = [4/7,4/7,13];
fiber_mas = 42; % diameter 
n = 150;
fx = normrnd(0,4/7,n,1);
fy = normrnd(0,4/7,n,1);
fz = normrnd(0,13,n,1);
rndpos = [fx,fy,fz];% global position offset in microns (x,y,z)

dof = 0; %Polychromatic depth of focus in microns

%-----Stellar Inputs-----%
units = 'counts'; % what physical units do you want? 'energy or counts'
type = 'M0V'; %Spectral type
v_mag = 7;%9+1.848; %Magnitude 
%V-I --> [F5 = 0.506, G0 = 0.664, G5 = 0.738,K0 = 0.853,K5 = 1.246, M0 = 1.848, M4=2.831, M7=4.52];
vsini = 2.5;
epsilon = 1; 
rv = 0; 
scale =1; %etalon spectrum upscaling factor


%-----Observation Conditions-----%
z_deg = 45; % %Zenith angle use 10-45-60 as defaults
AO = 2; % AO System & Seeing %1 = FLAO, 2 = SOUL 1'' , 3 =SOUL 0.8'', 4= SOUL 0.6'
Tellurics =true; % true or false 
SkyBack = false;

%-----Misc-----%
FSR = 10; %Etalon frequency GHz 
EtalonType = 'science'; %"science" or "ref" or "fullref" or "gauss" or delta

% starinputs = [] 
% etaloninputs = [a b c d]
% classinputs

% ---------------------------------------------------------------%
tic
% mag = [7:0.5:14]; %M0V
% ADC = [0:5:60];
% deg = [0:5:60];
output = [];

for ii = 1%1:length(mag)
% fiberpos = rndpos(ii,:);
% adc = ADC(ii); 
% z_deg = deg(ii);
% v_mag = mag(ii)+0.664;
%polarization = [1,ii,1];
dname = 'Science';
fname = ['TraceTestEntWindow']; %sets name
SimulationMain
% output = [output out];
end
toc

%% ------------------Change log -------------------------------%     
% changed the upsample factor and inverted the polynomial
% fixed the speed of light in Etalon
% set optical and pixel sampling to maximum 
% 8/22/17
% no upscaling of stellar spectrum for now
% no 2D wavelength solution yet
% generate etalon at appropriate sampling for convolution
% 10 Ghz 
% 8/23/17 
% New copy of code
% implementing 2D solution with interpolant
% 9/20/17
% New copy of code from V8
% updating all classes to include constructors
% introducing new class called fibercoupling which handles air-fiber
% interface
% 10/31/17
% V10 has fathom throughput, upsampled throughput, and real ADC curves