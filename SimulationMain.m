%RVMainTest

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
global adcfile

%% Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end

colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};
colors{11} =[0 0.3 0];
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
zmx_detector = [pathprefix '/RefFiles/Zemax/v12_22footprint.txt'];
gratingfile = [pathprefix '/RefFiles/Curves/GratingEff_new.mat'];
calibfile = [pathprefix '/Calibration/Lorentzian.mat'];
zernikefile = [pathprefix '/RefFiles/AO/zernike_index.mat'];
adcfile = [pathprefix '/RefFiles/ADC/ADC_coeff_radial.mat'];
%end of path block

%% ------------Initialise Objects------------ %%
%% Spectra Objects - stellar, etalon, flats. Creates instance of Simulation object called Target
Target = Simulation(type,v_mag,epsilon,vsini,TestRV,z_deg,spectrum_units,fiber_mas);

%Spectrum values are overwitten from this point on by Tellurics, Instrument profile etc.
if Tellurics == true 
Target.Tellurics(Target.Tellurics < 0)=1e-30; %tellurics are interpolated below 0% reset min value to 0    
Target.Spectrum(:,1) = Target.Tellurics(:,1).*Target.Spectrum(:,1); %multiply by tellurics (comment out to remove tellurics)
Target.Counts(:,1) = Target.Tellurics(:,1).*Target.Counts(:,1); %multiply by tellurics (comment out to remove tellurics)
end


if SkyBack == true

    [Resamp_Sky] = Target.Resamp_Extrap(Target.SkyBackground(:,1),Target.SkyBackground(:,2),Target.Wavelength); %Interpolate and extrapolate Energy and Wavelength onto DsWavelength;
    %temp = Target.Counts; %photons/sec/micron/m^2 for Skybackground figure
    Target.Counts = Target.Counts + Resamp_Sky;

end

Target = Target.VacShift; %Vacuum wavlength shift DsWavelength and Wavlength

%Multiply Spectrum by the collecting area of the LBT
Target.Counts(:,1) = Simulation.LBT_CollectingArea(Target.Counts(:,1));
%Multiply the delta lambda from counts to remove per wavelength unit. (technically it was counts/micron)
Target.Counts(:,2) = Target.Counts(:,1) * mean((Target.DsWavelength(2:end)-Target.DsWavelength(1:end-1)));

%Photon Noise requirements Butler and Marcy
% [out] = SignaltoNoise(Target.DsWavelength,Target.Counts(:,2));
% return
%Calibration Spectra
CalibSpec = Etalon(EtalonType,FSR,scale,TestRV);
% [out] = SignaltoNoise(CalibSpec.Wavelength/10000,CalibSpec.Counts(:,1));
% return

FlatSpec = Flat(scale);
% End of the spectra block
%% Fiber Object
SMF = FiberCoupling(wfe,adc,fiberpos,dof);

%% WFS Object 
LBTI_AO = WFSensor(AO,entwindow); % initialize WFS object  

%Multiply the WFS channel throughput by the spectrum
[Target.AOSpectrum(:,2),Target.AOSpectrum(:,1)] = Simulation.IncludeThroughput...
(LBTI_AO.FinalTput(:,1)/1000,Target.DsWavelength,LBTI_AO.FinalTput(:,2),Target.Counts(:,1));
Target.AOSpectrum(:,2) = Simulation.LBT_CollectingArea(Target.AOSpectrum(:,2)); %multiply counts by collecting are of the LBT
Target.IntegratedCounts = trapz(Target.AOSpectrum(:,1),Target.AOSpectrum(:,2)); %integrate counts using trapz method 
Target = CalculateStrehl(Target,AO); %calculate the Strehl ratio based on the total integrated counts
[Target.StrehlRatio(:,3)] = Simulation.rescale_SR(Target.StrehlRatio(:,2),Target.Airmass(1,2)); %Rescale Strehl Ratio according to airmass

% end of WFS channel block
%% iLocater object
FWR = iLocater(entwindow,band,'FWR',SMF.Rho,SMF.Wavelength,[0,0.5,0]);
FWR = SetInstrumentSR(Target.StrehlRatio,FWR);
FWR = Optical_Path(FWR); %Set Optical Path parameters
FWR = Path_Multiply(FWR); %Mulitply curves in optical path
% FWR = CustomizeBandpass(FWR);%Customize Bandpass (all have defaults preloaded)
FWR = Trim_Throughput(FWR);%Trim to custom Bandpass (1,:)
FWR = Integrated_Transmission(FWR); %Calc Integrated Throughput (no R6 grating)

ilocater = iLocater(entwindow,band,StarlightPath,SMF.Rho,SMF.Wavelength,polarization); %Name the Throughput object (no constructor right now)
ilocater = SetInstrumentSR(Target.StrehlRatio,ilocater); %Grab the Strehl ratio from the Simulation object and multiply into FIBERLINK!
ilocater = Optical_Path(ilocater); %Set Optical Path parameters
ilocater = Path_Multiply(ilocater); %Mulitply curves in optical path
ilocater = CustomizeBandpass(ilocater);%Customize Bandpass (all have defaults preloaded)
ilocater = Trim_Throughput(ilocater);%Trim to custom Bandpass (1,:)
ilocater = Integrated_Transmission(ilocater); %Calc Integrated Throughput (no R6 grating)

% out(:,1) = SMF.Wavelength;
% out(:,2) = SMF.Rho;
return

%% Calibration object
%Initialise dedicated object for calibration fiber
CalSource = Calibration(band);
CalSource = CustomizeBandpass (CalSource,[900 1350]);%Customize Bandpass (all have defaults preloaded)
CalSource = Trim_Throughput(CalSource);%Trim to custom Bandpass (1,:)
CalSource = Integrated_Transmission(CalSource);
% End of Instrument throughput block

%% ------------Combine instruments with spectra------------ %%
%% 1D spectrum (for acquisition Camera usage)

%Combine the Simple instrument with the spectrum (no R6 Grating)
[Target.InstrumentSpectrum(:,2),Target.InstrumentSpectrum(:,1)] = Simulation.IncludeThroughput...
    (ilocater.FinalTput(:,1)/1000,Target.DsWavelength,ilocater.FinalTput(:,2),Target.Counts(:,1));
Target.Bandpass = ilocater.Bandpass/1000;

%Combine the perfect instrument (100% throughtput) with the spectrum (no R6 Grating)
[Target.InstrumentSpectrum(:,3),~] = Simulation.IncludeThroughput(ilocater.FinalTput(:,1)/1000,...
    Target.DsWavelength,ones(length(ilocater.FinalTput(:,2))),Target.Counts(:,1));

%Combine the femptowatt reciever with the stellar spectrum.
[Target.FWRSpectrum(:,2),Target.FWRSpectrum(:,1)] = Simulation.IncludeThroughput...
    (FWR.FinalTput(:,1)/1000,Target.DsWavelength,FWR.FinalTput(:,2),Target.Counts(:,1));

%Calculate energy in Watts for each of the above
[Energy(:,1)] = Spectra.Counts2Energy (Target.InstrumentSpectrum(:,1),Target.InstrumentSpectrum(:,2));
[Energy(:,2)] = Spectra.Counts2Energy (Target.InstrumentSpectrum(:,1),Target.InstrumentSpectrum(:,3));
[EnergyFWR(:,1)] = Spectra.Counts2Energy (Target.FWRSpectrum(:,1),Target.FWRSpectrum(:,2));

Target.IntegratedEnergy(2,1) = trapz(Target.InstrumentSpectrum(:,1),Energy(:,2));
Target.IntegratedEnergy(3,1) = trapz(Target.InstrumentSpectrum(:,1),Energy(:,1));
Target.IntegratedEnergy(4,1) = trapz(Target.FWRSpectrum(:,1),EnergyFWR(:,1));

Target.IntegratedCounts(2,1) = trapz(Target.InstrumentSpectrum(:,1),Target.InstrumentSpectrum(:,3));
Target.IntegratedCounts(3,1) = trapz(Target.InstrumentSpectrum(:,1),Target.InstrumentSpectrum(:,2));
Target.IntegratedCounts(4,1) = trapz(Target.FWRSpectrum(:,1),Target.FWRSpectrum(:,2));

if strcmp(ilocater.PathName,'Spectrograph') == 0 %only Spectrograph will proceed from here
    return
end
%% Break the spectrum out into orders according to the throughput. 

%Prepare Spectrum to broken into orders by combining the E2E insturment with grating (R6 Grating with all simulated orders)
ilocater = Include_Grating(ilocater); %Include grating into throughput object by creating a new property Instrument.Orders
CalSource = Include_Grating(CalSource);%Include grating into throughput object by creating a new property Calibration.Orders

if Throughput == 1
    StarEfficiency = ilocater.Orders{1,1};
    CalEfficiency = CalSource.Orders{1,1};
elseif Throughput == 0
    StarEfficiency = ones(size(ilocater.Orders{1}));
    CalEfficiency  = ones(size(ilocater.Orders{1}));
end


for ii = 1:size(ilocater.Orders{2},2)-3 %(obj.Orders is now in cell array format due to grating)
        % Spectral Orders with full instrument profile
        [Target.SpectrumOrders{2}(:,ii),Target.SpectrumOrders{1}(:,ii)] = Simulation.IncludeThroughput(ilocater.Orders{1,2}(:,ii)./1000,...
            Target.DsWavelength,StarEfficiency(:,ii),Target.Counts(:,2));
        
        % Etalon Orders with full instrument profile
        [Target.EtalonOrders{2}(:,ii),Target.EtalonOrders{1}(:,ii)] = Simulation.IncludeThroughput(CalSource.Orders{1,2}(:,ii)./1000,...
            CalibSpec.Wavelength./1e4,CalEfficiency(:,ii),CalibSpec.Counts(:,1));
        
        % Flat Orders with full instrument profile
        [Target.FlatOrders{2}(:,ii),Target.FlatOrders{1}(:,ii)] = Simulation.IncludeThroughput(CalSource.Orders{1,2}(:,ii)./1000,...
            FlatSpec.Wavelength,CalEfficiency(:,ii),FlatSpec.Counts(:,1));
        
       % Reference Orders with 100% throughput
         [NumTput{2}(:,ii),NumTput{1}(:,ii)] = Simulation.IncludeThroughput(ilocater.Orders{1,2}(:,ii)./1000,...
         Target.DsWavelength,ones(size(ilocater.Orders{1})),Target.Counts(:,2));
end

% All orders can be found in Simulation object: SpectrumOrders,EtalonOrders,FlatOrders

            
%Now there is a fully prepared Simulation Object.
%% Scale the calibration Spectra in relation to the stellar continuum
%[Target.EtalonOrders{2}]=Etalon.Scale2Star(Target.EtalonOrders{2},Target.SpectrumOrders{2});
%[Target.FlatOrders{2}]=Flat.Scale2Star(Target.FlatOrders{2},Target.SpectrumOrders{2});
%Target.EtalonOrders{2}= Target.EtalonOrders{2}./4;

if Throughput == 1
WLS =0.1e-3*1000; %0.1mW/nm converted to mW/um (multiply by 1000)
WLS = WLS/1.2e-19; %ph/s/um;
WLS = WLS *1.2e-6;% 2.55e-6 photons/s/bin (assuming same sampling as a star)
CT = 0.5^10*0.1*0.01*0.01;%calibration transmission (through fiber chains)
[Target.EtalonOrders{2}]=(Target.EtalonOrders{2}*CT*WLS); %WLS*transmission through spectrograph 
elseif Throughput == 0
% do nothing
end

%% ------------Put down spectra on detector------------ %%
%% Minimize Memory usage before parallel loop
% Create wavelength solution or load it from previous solution
if wavesolution == 1
    
    if cheby ==1
        for ii = 1:3
            [p1{ii},ret{ii},chebs{ii}] = WaveSolution(zmx_detector,ii);
%             save chebycoeffs2 chebs p1 ret
        end
    end
    [~,order_coeff,wave_coeff]=Simulation.ConvertZemaxFormat(zmx_detector);
     save polycoeffs2 order_coeff wave_coeff
else
    
    load chebycoeffs2
    load polycoeffs2
end

Orders_Name = {'Etalon','Science', 'Flat'};

for jj = 1:3 %number of order variables (Etalon , Flat , Star)
    
    if strcmp(Orders_Name{jj},'Etalon')==1
        inputOrders = Target.EtalonOrders{2};
        inputWavelength = Target.EtalonOrders{1};
    elseif strcmp(Orders_Name{jj},'Science')==1
        inputOrders = Target.SpectrumOrders{2};
        inputWavelength = Target.SpectrumOrders{1};
    elseif strcmp(Orders_Name{jj},'Flat')==1
        inputOrders = Target.FlatOrders{2};
        inputWavelength = Target.FlatOrders{1};
    end
    
    
    [outputWavelength, outputOrders] = Simulation.trimOrders(inputWavelength, inputOrders, wave_coeff);
   
    if strcmp(Orders_Name{jj},'Science')==1
        [Num_outputWavelength, Num_outputOrders] = Simulation.trimOrders(NumTput{1},NumTput{2}, wave_coeff);
        
        %% Numerical Throughput
        for ii = 1:size(outputOrders,2)
            NumTput{3}(:,ii) = outputOrders(:,ii)./Num_outputOrders(:,ii);
        end
        NumThroughput = mean(mean(NumTput{3}));
    end
    
    % if trim orders fails check to make sure zemax file orders are in
    % correct order (not middle out) this is a manual fix
    if strcmp(Orders_Name{jj},'Etalon')==1
        Target.EtalonOrders{1} = outputWavelength;
        Target.EtalonOrders{2} = outputOrders;
    elseif strcmp(Orders_Name{jj},'Science')==1
        Target.SpectrumOrders{1} = outputWavelength;
        Target.SpectrumOrders{2} = outputOrders;
    elseif strcmp(Orders_Name{jj},'Flat')==1
        Target.FlatOrders{1} = outputWavelength;
        Target.FlatOrders{2} = outputOrders;
    end
end


% central trace is indexed to jj == 1, the others are unkown currently. i
% think jj == 3 is the "bottom".

% Reference values for the fits header
ind = find(Target.StrehlRatio(:,1)>= 1.12,1,'First');
mid_SR = Target.StrehlRatio(ind,3);
Rho = mean(ilocater.FiberLink(:,4));

if mode == 0  % Science frame
    
    headerinfo = {'SpType', Target.SpType, ' ';...
        'Vmag', Target.Vmag,' ';...
        'Rmag', Target.Rmag,' ';...
        'Imag', Target.Imag,' ';...
        'RV', Target.RV,'Injected RV (m/s)';...
        'Vsini', Target.Vsini,' ';...
        'Tlrcs', Tellurics ,'Are tellurics included 1/0';...
        'Backgnd', SkyBack,'Is skybackground included? 1/0';...
        'Airmass', Target.Airmass(1,3),' ';...
        'Strehl', mid_SR,'Average Strehl';...
        'Rho', Rho,'Average SMF coupling';...
        'xSamp', scale,'Upsampled scale factor';...
        'CalTrace','Etalon','Parameters Below';...
        'FSR',CalibSpec.FSR,'GHz';...
        'Finesse',CalibSpec.Finesse,' ';...
        'AO',LBTI_AO.AO,'1,2,3 = SOUL, 0 = LBTIAO';...
        'SMFPosX',fiberpos(1),' ';...
        'SMFPosY',fiberpos(2),' ';...
        'SMFPosZ',fiberpos(3),' ';...
        'ADCPos',adc,'max dispersion';...
        'Tput',NumThroughput,'mean Throughput';...
        'DOP',polarization(1,1),'Degree of polarization';...
        'Pfrac',polarization(1,2),'Fraction of P-pol';...
        };
    
    fitsname = [pathprefix,'/Output/Science/',dname,fname,'.fits'];
    trace = {'Science','Science','Etalon'}; %(Bottom, Top, Middle)
    
elseif mode == 1
    
    headerinfo = {'Type', 'Flat', 'Spectral type';...
        'xSamp', scale,'Upsampled scale factor';...
        };
    
    fitsname = [pathprefix,'/Output/Flat/',dname,fname,'.fits'];
    trace = {'Flat','Flat','Flat',}; %(Bottom, Top, Middle)
     
elseif mode == 2
    headerinfo = {'Type', 'Gauss', 'Spectral type';...
        'xSamp', scale,'Upsampled scale factor';...
        'RV', CalibSpec.RV,'Injected RV (m/s)';...
        };
    
    fitsname = [pathprefix,'/Output/Etalon/',dname,fname,'.fits'];
    trace = {'Etalon','Etalon','Etalon',}; %(Bottom, Top, Middle)
    
elseif mode == 3
    headerinfo = {'Type', 'Testing', 'Spectral type';...
        'xSamp', scale,'Upsampled scale factor';...
        'RV', Target.RV,'Injected RV (m/s)';...
        };
    
    fitsname = [pathprefix,'/Output/Testing/',dname,fname,'.fits'];
    trace = {'Science','Etalon','Flat',}; %(Bottom, Top, Middle)
    
end
toc

for jj = tracenum
    if strcmp(trace{jj},'Etalon')==1
        inputspec = Target.EtalonOrders{2};
        inputwave = Target.EtalonOrders{1};
        isStar = 0;
        
    elseif strcmp(trace{jj},'Science')==1
        inputspec = Target.SpectrumOrders{2};
        inputwave = Target.SpectrumOrders{1};
        isStar = 1;
    elseif strcmp(trace{jj},'Flat')==1
        inputspec = Target.FlatOrders{2};
        inputwave = Target.FlatOrders{1};
        isStar = 0;
    end
    fprintf('\nStarting trace %i using %s spectrum \n',jj, trace{jj})
    if parflag == true
        parfor ii = 1:norders
            
            fprintf('Computing order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(inputwave(:,ii),inputspec(:,ii),wave_coeff(ii,:,jj),scale);
            Detector(:,:,ii,jj) = Simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))
            
        end
    else
        for ii = 1:norders
            
            fprintf('Computing order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(inputwave(:,ii),inputspec(:,ii),wave_coeff(ii,:,jj),scale);
            Detector(:,:,ii,jj) = Simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))

        end
    end
    
    
end
disp('Writing to .fits file')
% Target.DetectorFace = squeeze(sum(Detector(:,:,:,jj),3)); % for just the jth orderline 
Target.DetectorFace = sum(sum(Detector,4),3);
% Target.DetectorFace(:,:,4) = sum(Target.DetectorFace,3);
% Target.DetectorFace = fliplr(rot90(Target.DetectorFace,2));
% TempDetectorFace = Target.DetectorFace(:,:,jj);

Simulation.WriteFits(fitsname,Target.DetectorFace,headerinfo)

save backupArray Target

disp('done')
toc






