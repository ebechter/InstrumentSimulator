% Simulation Main Script
% clear up the workspace
clear; close all; clc
addpath(genpath(pwd))


% First thing is to check user inputs from the wrapping script.

% Create initial instances of classes (objects)

% mode 2 (SES) part of simulation






% 
% for ii = tracenum
%         
%         
%         trace{1} = curve{1};
%         trace{2} = curve{2};
%         trace{3} = curve{1};
%         
%     end
    
    
    
    
    
    
    
    
    
    
    
    
%     
% end
% 
% 
% Spectroscopy Modes
% 
% I want mode x and traces (default all 3)  [ 1 -> 3]
% 
% mode 1
% 
% full star -> detector
% 
% mode 2
% 
% etalon -> calibration throughput -> spectrograph
% 
% mode 3
% 
% flat
% 
% 
% 
% 
% 
% 
% Operation for traces [1 -> 3]
% 
% trace{1} = curve 1: typical star thing
% 
% trace{2} = curve 2: typical etalon
% 
% trace{3} = curve 1: repeat star
% 
% 
% 
% 
% curve 1:
% 
% atmosphere? yes/no
% throughput? yes/no
% 
% source: star, atmosphere: yes ; throughputs: telescope, lbti, acquisition_camera, spectrograph
% 
% 
% 
% 
% 
% 
% 
% mode 2:
% 
% sources: star, etalon, star
% 
% atmosphere: yes, no, yes
% 
% throughputs: telescope, lbti, acquisition_camera, spectrograph ; calibration; telescope, lbti, acquisition_camera, spectrograph
% 
% mode 3:
% 
% sources: etalon, etalon, etalon
% 
% atmosphere: no, no, no
% 
% throughputs: calibration, calibration, calibration
% 
% 
% mode 4:
% 
% sources: flat, flat, flat
% 
% atmosphere: no, no, no
% 
% throughputs: calibration, calibration, calibration



curve{1}.source = 'Star';
curve{1}.atmosphere = 1;
curve{1}.throughput = {'lbt','lbti','fiber','spectrograph'};

curve{2}.source = 'Etalon';
curve{2}.atmosphere = 0;
curve{2}.throughput = {'calibration','spectrograph'};

curve{3} = curve{1};



%========== Instanciate objects ===========%
% tracenum = [1,2,3];
tracenum = 1;
for ii = tracenum
    
    %========== Source Options ===========%
    
    if strcmp('Etalon', curve{ii}.source) == 1 && exist('etalon','var') == 0
        % make an etalon
        etalon = Etalon();
        
    elseif strcmp('Star', curve{ii}.source) == 1 && exist('star','var') == 0
        % make a star
        star = Star();
        
    elseif strcmp('Flat', curve{ii}.source) == 1 && exist('flat','var') == 0
        % make a flat spectrum
        flat = Flat();
    end
    
    %========== Atmosphere Options ===========%
    
    if curve{ii}.atmosphere == 1 && exist('atmosphere','var') == 0
        atmosphere = Atmosphere();
    end
    
    %========== Throughput Options ===========%
    components =[];
    
    if nonzeros(strcmp('spectrograph', curve{ii}.throughput)) == 1 && exist('spectrograph','var') == 0
        % make the spectrograph throughput
        spectrograph = Spectrograph();
        
    end
    
    if nonzeros(strcmp('lbt', curve{ii}.throughput)) == 1 && exist('lbt','var') == 0
        % make the lbt throughput
        lbt = Imager('LBT');
        components = [components, lbt];
        
    end
    
    if nonzeros(strcmp('lbti', curve{ii}.throughput)) == 1 && exist('lbti','var') == 0
        % make the l throughput
        lbti = Imager('LBTI');
        components = [components, lbti];
    end
    
%     if nonzeros(strcmp('lbti', curve{ii}.throughput)) == 1 && exist('lbti','var') == 0
%         % make the l throughput
%         lbti = Imager('LBTI');
%         components = [components, lbti];
%     end
    
    if nonzeros(strcmp('fiber', curve{ii}.throughput)) == 1 && exist('fiber','var') == 0
        % make the l throughput
        fiber = Imager('FIBER');
        components = [components, fiber];
    end
end

% Combine imager throughputs  
%  [finalThroughput] = combinedImagerThroughput(objects);
 comTput{1} = components(1).finalThroughput;
 
 for ii = 2:size(components,2) %loop over cellarry size (i.e. surfaces)
     [comTput{ii}(:,1),comTput{ii}(:,2)] = Instrument.multiply_curves...
         (comTput{ii-1}(:,1),comTput{ii-1}(:,2),...
          components(ii).finalThroughput(:,1),components(ii).finalThroughput(:,2));    
 end




if exist('atmosphere','var')
%     [addedatmosphere] = Simulation.AddAtmosphere(star,atmosphere);
[skygrid] = Simulation.resampleToGrid(atmosphere.skyback(:,1),atmosphere.skyback(:,2),star.dsWavelength);

atmosphere.skyback = [star.dsWavelength,skygrid];

addsky = star.counts + skygrid;

[tell_grid] = Simulation.resampleToGrid(atmosphere.telluric(:,1),atmosphere.telluric(:,2),star.dsWavelength);

addtell = addsky.*tell_grid;
else
    
    addedatmosphere = star;
end


[throughput_grid] = Simulation.resampleToGrid(comTput{end}(:,1)*1e-3,comTput{end}(:,2),star.dsWavelength);

addthroughput = addtell.*throughput_grid;


% 
% [modifiedspectrum] = Simulation.MultiplyThroughput(addedatmosphere, combinedthroughput);
% 
% if exist('spectrgraph','var')
%     [modifiedspectrum] = Simulation.MultiplyThroughput(addedatmosphere, combinedthroughput);
% end
% 


% Pass them to the Simulation class

% simulation class








