function SimulationMainP(parallelInfo,runInfo, specInfo,fitsname, ... 
                    wfe,starInfo,source,polarization,SpecOrImager,conditions, headerinfo,persistentSource)

poolobj = gcp('nocreate');
nOrders = specInfo{2};
cheby = specInfo{3};
p1 = specInfo{4};
ret = specInfo{5};
tracenum = specInfo{6};
order_coeff = specInfo{7};
wave_coeff = specInfo{8};
zenith = conditions{1};
seeing = conditions{2};
entWindow = conditions{3};
aoType = conditions{4};
scale = runInfo{2};
persistence = persistentSource.persistence; 




if parallelInfo{2}
    usep = 'yes';
    if isempty(poolobj)
        %no pool exists and I want one
        parpool(parallelInfo{1});
        
    else
        %just run it b/c a pool already exists
    end
else
    usep = 'no';
    delete(poolobj); % should shut down the current pool or does nothing if there is no pool
end

fprintf('-----\n%s simulation\nversion: %s\nfootprint: %s\nupscale factor: %i\n-----\nparallel settings\nParallel: %s\nnum. of cores: %i\n-----\n'...
        ,SpecOrImager, runInfo{1}, specInfo{1}, runInfo{2}, usep, parallelInfo{1})

for ii = tracenum
    fprintf('creating simulation objects...\n')
    %========== Source Options ===========%
    
    if strcmp('etalon', source{ii}.name) == 1 && exist('etalon','var') == 0
        % make an etalon
        etalon = Etalon(scale);
        
    elseif strcmp('star', source{ii}.name) == 1 && exist('star','var') == 0
        % make a star

        star = Star(starInfo{1},starInfo{2},starInfo{3},starInfo{4},starInfo{5},starInfo{6});
        
    elseif strcmp('flat', source{ii}.name) == 1 && exist('flat','var') == 0
        % make a flat spectrum
        flat = Flat(scale);
    
    elseif strcmp('superk', source{ii}.name) == 1 && exist('superk','var') == 0
        superk = SuperK();
    end
    
    
    %========== Persistent Source Options ===========%
 
    if strcmp('etalon', persistentSource.name) == 1 && exist('persistentEtalon','var') == 0
        % make an etalon
        persistentEtalon = Etalon(3,persistentSource.rv); % scale of 3 is fine here - interpolated down automatically later.
        
    elseif strcmp('star', persistentSource.name) == 1 && exist('persistentStar','var') == 0
        % make a star

        persistentStar = Star(persistentSource.spType,persistentSource.vmag,persistentSource.epsilon,...
            persistentSource.vsini,persistentSource.rv,persistentSource.units);
        
    elseif strcmp('flat', persistentSource.name) == 1 && exist('persistentFlat','var') == 0
        % make a flat spectrum
        persistentFlat = Flat(scale);
    
    elseif strcmp('superk', persistentSource.name) == 1 && exist('persistentSuperk','var') == 0
        persistentSuperk = SuperK();
    end
    
    
    
    %========== Atmosphere Options ===========%
    
    if source{ii}.atmosphere == 1 && exist('atmosphere','var') == 0
        atmosphere = Atmosphere();
    end
    
    %========== Throughput Options ===========%
    if exist('imagerComponents','var')==0 && exist('imagerComponents','var') ==0
        imagerComponents =[];
    end
    
    if exist('AO_list','var')==0 && exist('aoComponents','var') ==0
        aoComponents = [];
    end
    
    if any(strcmp('spectrograph', source{ii}.throughput)) == 1 && exist('spectrograph','var') == 0
        % make the spectrograph throughput
        spectrograph = Spectrograph(polarization);
    end
    
    if any(strcmp('lbt', source{ii}.throughput)) == 1 && exist('lbt','var') == 0
        % make the lbt throughput
        lbt = Imager('LBT');
        imagerComponents = [imagerComponents, lbt];
    end
    
    if any(strcmp('lbti', source{ii}.throughput)) == 1 && exist('lbti','var') == 0
        % make the l throughput
        lbti = Imager('LBTI');
        imagerComponents = [imagerComponents, lbti];
    end
    
    if any(strcmp('fiberCh', source{ii}.throughput)) == 1 && exist('fiberCh','var') == 0
        % make the l throughput
        fiber = Imager('FIBER');
        imagerComponents = [imagerComponents, fiber];
    end
    
    if any(strcmp('imageCh', source{ii}.throughput)) == 1 && exist('fiberCh','var') == 0
        % make the l throughput
        andor = Imager('ANDOR');
        imagerComponents = [imagerComponents, andor];
    end
    
    if any(strcmp('quadCh', source{ii}.throughput)) == 1 && exist('quadCh','var') == 0
        % make the l throughput
        quad = Imager('QUADCELL');
        imagerComponents = [imagerComponents, quad];
    end
    
    if any(strcmp('wfc', source{ii}.throughput)) == 1 && exist('wfc','var') == 0
        % make the l throughput
        wfc = Imager('WFC');
        imagerComponents = [imagerComponents, wfc];
    end

    if source{ii}.AO == 1 && exist('lbti_ao','var') == 0
        % make the l throughput
        lbti_ao = AO([aoType entWindow]);
        aoComponents = [aoComponents, lbti_ao];
    end
    
    if any(strcmp('filter', source{ii}.throughput)) == 1 && exist('filter','var') == 0
        % make the lbt throughput
        filter = Imager('Filter');
        imagerComponents = [imagerComponents, filter];
    end
    
    if any(strcmp('flat', source{ii}.throughput)) == 1 && exist('flat','var') == 0
        % make the lbt throughput
        flat = Imager('flat');
        imagerComponents = [imagerComponents, flat];
    end
    
end

simulation = Simulation(scale); % upscaling simulation
spectral_cell = cell(3,1);


% this decides how to combine throughput terms with sources and atmospheres for each of the traces.

for ii = tracenum

    if strcmp(source{ii}.name,'star')
        
        if strcmp(SpecOrImager,'Spectrograph') == 1
            spectral_cell{ii} = Simulation.addStar(spectrograph.maxR,spectrograph.pixSamp, simulation.scale, ...
                star.wavelength,star.spectrum,star.rv,0); % spectrum units in W/m^2/um
            
            
                
                if exist('persistentStar','var') == 1
                    
            pStar = Simulation.addStar(spectrograph.maxR,spectrograph.pixSamp, simulation.scale, ...
                persistentStar.wavelength,persistentStar.spectrum,persistentStar.rv,persistence); % spectrum units in W/m^2/um
            spectral_cell{ii} = Simulation.CombineSpectra(spectral_cell{ii},pStar);
            clear pStar
                end
                if exist('persistentEtalon','var') == 1
                    
                    pEtalon(:,1) = persistentEtalon.wavelength;
                    pEtalon(:,2) = persistentEtalon.counts* max(spectral_cell{ii}(:,2)) * persistence; % for now need to scale etalon to star (roughly)

                    spectral_cell{ii} = Simulation.CombineSpectra(spectral_cell{ii},pEtalon);
                    clear pEtalon persistentEtalon
                end
        elseif strcmp(SpecOrImager,'Imager') == 1
            
                spectral_cell{ii}(:,1) = star.wavelength;
                spectral_cell{ii}(:,2) = Star.energy2Counts(star.wavelength,star.spectrum);
        end
        
        if source{ii}.atmosphere == 1
            
            spectral_cell{ii}(:,2) = Simulation.addAtmosphere(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2), ...
                atmosphere.telluric, atmosphere.skyback);
            
        end
        
        spectral_cell{ii}(:,2) = Simulation.addCollectingArea(spectral_cell{ii}(:,2),lbt.apDiameter,lbt.blockFrac); % spectrum units in....
        spectral_cell{ii}(:,2) = Star.fluxDenToflux(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2)); % spectrum units in....
        
        if source{ii}.AO == 1
            fprintf('calculating strehl ratio...')

            % if we want ao, first check if its already been made
            if exist('AO_throughput','var') == 0
                
                % if not, make it
                [AO_throughput] = Simulation.combineImagerThroughput(aoComponents);
            end
            
            % trim source to WFS band (could change between traces)
            [wfsWave,wfsSpec] = Simulation.trimToBand(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2),lbti_ao.bandPass*1e-3);
            
            % compute counts on wfs.
            AOFlux = Simulation.resampleToGrid(AO_throughput(:,1)*1e-3,AO_throughput(:,2),wfsWave);
            wfsCounts = sum(wfsSpec.*AOFlux);
            
            strehlR = Simulation.calculateStrehlRatio(aoType,seeing,wfsCounts,zenith);
            fprintf('%s \n',char(hex2dec('2713')))

        end
        
        if isempty(source{ii}.throughput) == 0
            
            if exist('starThroughput','var') ==0
                
                [starThroughput,tputProg] = Simulation.combineImagerThroughput(imagerComponents);
            end
            
            throughputGrid = Simulation.resampleToGrid(starThroughput(:,1)*1e-3,starThroughput(:,2),spectral_cell{ii}(:,1));
            
            if source{ii}.tput_flag == 1
            spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2).*throughputGrid;
            else 
            spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2);
            end

        end
        
        
        if any(strcmp('fiberLink', source{ii}.throughput)) == 1
            
            fprintf('calculating fiber coupling...')
            
            strehl = strehlR; 
            strehlR(:,1) = strehlR(:,1)*1000;% need to work in nm for fiber coupling
            wfef = [0,0,0,0,0,0,0,0]; % complicated input (user can specify wave front error)
            adc = 0; % zenith angle for ads 0-60 in steps of 5
            fiberpos = [0,0,0]; % global position offset in microns (x,y,z)
            dof = 0; % depth of focus (not sure if used yet)
            bandPass = 965:1305; % typical bandpass for fiber coupling
            smfLink = FiberLink(wfef,adc,fiberpos,dof,strehlR,bandPass);
            
            % indlcue all fiber coupling losses (Strehl ratio, mismatch,
            % adc, fiber offsets etc) specLink is for spectrograph path
            
            %multiply the star by fiber throughput
            throughputGrid = Simulation.resampleToGrid(smfLink.specLink(:,1)*1e-3,smfLink.specLink(:,2),spectral_cell{ii}(:,1));
            spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2).*throughputGrid;
            
            %record throughput progression
            new_y = Simulation.resampleToGrid(smfLink.specLink(:,1),smfLink.specLink(:,2),tputProg{1,end}(:,1));
            tempTput(:,1) = tputProg{1,end}(:,1);
            tempTput(:,2) = new_y .* tputProg{1,end}(:,2);
            
            %[tempTput] = Simulation.addSpecThroughput(,spectrograph.finalThroughput,nOrders);
            smfCell{1} = tempTput;
            smfCell{2,1} = smfLink.name;
            tputProg = [tputProg smfCell];
            clear specCell tempTput new_y
            fprintf('%s \n',char(hex2dec('2713')))

        end

    elseif strcmp(source{ii}.name,'etalon')

        spectral_cell{ii} = [etalon.wavelength etalon.counts];
        tputProg{1,1} = [etalon.wavelength*1e3,ones(size(etalon.counts))];
        tputProg{2,1} = 'EtalonOnes';
    
    elseif strcmp(source{ii}.name,'flat')
        
        spectral_cell{ii} = [flat.wavelength flat.counts];
        tputProg{1,1} = [flat.wavelength*1e3,ones(size(flat.counts))];
        tputProg{2,1} = 'FlatOnes';
        
    elseif strcmp('superk', source{ii}.name) == 1
        
        spectral_cell{ii}(:,1) = superk.wavelength;
        spectral_cell{ii}(:,2) = superk.counts;
        spectral_cell{ii}(:,2) = Star.fluxDenToflux(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2));
        
        if isempty(source{ii}.throughput) == 0
            
            if exist('starThroughput','var') ==0
                [Throughput,tputProg] = Simulation.combineImagerThroughput(imagerComponents);
                throughputGrid = Simulation.resampleToGrid(Throughput(:,1)*1e-3,Throughput(:,2),spectral_cell{ii}(:,1));
                spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2).*throughputGrid;
                
            end
            
        end
        
    end
   
    if strcmp(SpecOrImager,'Imager')
        
        % trim wavelength band to instrument band or user input
        bandPass = imagerComponents(1,end).bandPass./1000; % hard cutoff for bandpass in microns
        [trimWave,trimCounts] = Simulation.trimToBand(spectral_cell{1}(:,1),spectral_cell{1}(:,2),bandPass);
        
        %calculate energy and counts
        simulation.totalCounts = sum(trimCounts);
        simulation.totalEnergy = Spectra.counts2Energy(trimWave,trimCounts);
        simulation.totalEnergy = sum(simulation.totalEnergy);
        
        %total = (max(Bandpass(1,:)*1000)-min(Bandpass(1,:))*1000);
        %IntTrans = trapz(tputProg{1,3}(:,1),tputProg{1,3}(:,2))/total;
        
        
        
        %Simulate a frame
        counts = simulation.totalCounts;
        pixelpitch = imagerComponents(1,end).pixelPitch;
        dimensions = imagerComponents(1,end).detectorDimensions;
        psf = imagerComponents(1,end).psf;
        FR = 1;
%         lbt.progressPlot
%         lbti.progressPlot
%         fiber.progressPlot
%         simulation.ProgressionPlot(tputProg)
        [frame] = Simulation.simulateImager(psf,dimensions,pixelpitch,counts/FR); % currently uses gaussian beam but should use FT code for PSF generation.
        fprintf('writing detector face to a fits file...')
        save ImageArray frame
        simulation.WriteFits(fitsname,frame,headerinfo)
        fprintf('%s \n',char(hex2dec('2713')))
        
        %Don't need to go any further for imagers.
        %clear nOrders ii order_coeff ret scale cheby chebs tracenum wave_coeff p1 parflag SpecOrImager total spectral_cell
        
    end
    
        % cross disperse spectrum into nOrders orders, trim down wavelengths
        
        spectral_cell{ii} = Simulation.Xdisperse(spectral_cell{ii},nOrders,wave_coeff);
        
        if any(strcmp('spectrograph', source{ii}.throughput)) == 1
            
            % include throughput of spectrograph 
            if (source{ii}.tput_flag) == 0
                
                for jj = 1:size(spectrograph.finalThroughput{1},2)
                    
                    spectrograph.finalThroughput{1}(:,jj) = ones(size(spectrograph.finalThroughput{1}(:,jj)));
                end
                spectral_cell{ii} = Simulation.addSpecThroughput(spectral_cell{ii},spectrograph.finalThroughput,nOrders);

            else
                
                spectral_cell{ii} = Simulation.addSpecThroughput(spectral_cell{ii},spectrograph.finalThroughput,nOrders);
            
            end
            
            for nn = 1:nOrders
                new_y = Simulation.resampleToGrid(spectrograph.finalThroughput{2}(:,nn),spectrograph.finalThroughput{1}(:,nn),tputProg{1,end}(:,1));
                tempTput{1}(:,nn) = tputProg{1,end}(:,1);
                tempTput{2}(:,nn) = new_y .* tputProg{1,end}(:,2);
            end
            
            %[tempTput] = Simulation.addSpecThroughput(,spectrograph.finalThroughput,nOrders);
            specCell{1} = tempTput;
            specCell{2,1} = spectrograph.name;
            tputProg = [tputProg specCell];
            clear specCell tempTput new_y tputProg
            
        end
end  

% At this point you have a spectral_cell of {trace}{order}(wavelength,counts)

% Now convolve each spectral order and clip on detector in parallel or serial.

fprintf('-----\nclipping spectrum to detector...\n')
for jj = tracenum
    
    fprintf('trace %i \n',jj)
    if parallelInfo{2} == true
        parallel_scale = simulation.scale;
        parallel_cell = spectral_cell{jj};
        
        parfor ii = 1:specInfo{2}
            
            fprintf('working on order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(parallel_cell{1}(:,ii),parallel_cell{2}(:,ii),wave_coeff(ii,:,jj),wfe,parallel_scale,ii);
            detector(:,:,ii,jj) = simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))
            
        end
        
        clear temp_cell
    else
        for ii = 1:specInfo{2}
            
            fprintf('working on order %i... ',ii)
            
            [OrderFlux{ii}, OrderWave{ii}] = Simulation.ConvolveOrder(spectral_cell{jj}{1}(:,ii),spectral_cell{jj}{2}(:,ii),wave_coeff(ii,:,jj),wfe,simulation.scale,ii);
            detector(:,:,ii,jj) = simulation.CliptoDetector(OrderFlux{ii}, OrderWave{ii},order_coeff(ii,:,jj),wave_coeff(ii,:,jj),cheby,p1{jj},ii);
            fprintf('%s \n',char(hex2dec('2713')))
            
        end
    end
    
    
end

fprintf('-----\n')
frame = sum(sum(detector,4),3);
clear detector

fprintf('writing detector face to a fits file...')
save BackupArray frame
simulation.WriteFits(fitsname,frame,headerinfo)
fprintf('%s \n',char(hex2dec('2713')))
end
