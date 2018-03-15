% simulation Main Script
% clear up the workspace
clear; close all; clc
addpath(genpath(pwd))
parflag = true;
scale = 1;
load polycoeffs2
load chebycoeffs2
nOrders = 2;
cheby=0;

curve{1}.source = 'star';
curve{1}.atmosphere = 1;
curve{1}.throughput = {'lbt','lbti','fiber','spectrograph'};
curve{1}.AO = 1;

curve{2}.source = 'etalon';
curve{2}.atmosphere = 0;
curve{2}.throughput = {'calibration','spectrograph'};
curve{2}.AO = 0;

curve{3} = curve{1};

%========== Instanciate objects ===========%
tracenum = 1;%[1,2,3];
for ii = tracenum
    
    %========== Source Options ===========%
    
    if strcmp('etalon', curve{ii}.source) == 1 && exist('etalon','var') == 0
        % make an etalon
        etalon = Etalon();
        
    elseif strcmp('star', curve{ii}.source) == 1 && exist('star','var') == 0
        % make a star
        star = Star();
        
    elseif strcmp('flat', curve{ii}.source) == 1 && exist('flat','var') == 0
        % make a flat spectrum
        flat = Flat();
    end
    
    %========== Atmosphere Options ===========%
    
    if curve{ii}.atmosphere == 1 && exist('atmosphere','var') == 0
        atmosphere = Atmosphere();
    end
    
    %========== Throughput Options ===========%
    if exist('star_components','var')==0
        star_components =[];
    end
    if exist('AO_list','var') ==0
        AO_list = [];
    end
    if any(strcmp('spectrograph', curve{ii}.throughput)) == 1 && exist('spectrograph','var') == 0
        % make the spectrograph throughput
        spectrograph = Spectrograph();
        
    end
    
    if any(strcmp('lbt', curve{ii}.throughput)) == 1 && exist('lbt','var') == 0
        % make the lbt throughput
        lbt = Imager('LBT');
        star_components = [star_components, lbt];
        
        if curve{ii}.AO == 1
            AO_list = [AO_list, lbt];
        end
    end
    
    if any(strcmp('lbti', curve{ii}.throughput)) == 1 && exist('lbti','var') == 0
        % make the l throughput
        lbti = Imager('LBTI');
        star_components = [star_components, lbti];
    end
    
    if any(strcmp('fiber', curve{ii}.throughput)) == 1 && exist('fiber','var') == 0
        % make the l throughput
        fiber = Imager('FIBER');
        star_components = [star_components, fiber];
    end
    
    if curve{ii}.AO == 1 && exist('lbti_ao','var') == 0
        % make the l throughput
        lbti_ao = Imager('LBTI_AO');
        AO_list = [AO_list, lbti_ao];
                    
    end
    
end

simulation = Simulation(scale); 
spectral_cell = cell(3,1);


% this decides how to combine throughput terms with sources and atmospheres for each of the traces. 
runDecisionTree()

% At this point you have a spectral_cell of {trace}{order}(wavelength,counts)

% Now convolve each spectral order and clip on detector in parallel or serial. 
for jj = tracenum

    fprintf('\nStarting trace %i \n',jj)
    if parflag == true
        parallel_scale = simulation.scale;
        parallel_cell = spectral_cell{jj};
        
        parfor ii = 1:nOrders
            
            fprintf('Computing order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(parallel_cell{1}(:,ii),parallel_cell{2}(:,ii),wave_coeff(ii,:,jj),parallel_scale);
            Detector(:,:,ii,jj) = Simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))
            
        end
        
        clear temp_cell
    else
        for ii = 1:nOrders
            
            fprintf('Computing order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(spectral_cell{jj}{1}(:,ii),spectral_cell{jj}{2}(:,ii),wave_coeff(ii,:,jj),simulation.scale);
            Detector(:,:,ii,jj) = Simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))

        end
    end
    
    
end
















% Fiber throughput and AO correction go into combine Imager curves.  

% Then include the spectrograph orders and throughput. 

% Derive wavelength solution map or load saved map. 

% Convolve and Clip loop. 











