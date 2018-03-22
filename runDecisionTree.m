for ii = tracenum
   
    if strcmp(curve{ii}.source,'star') 
        
        spectral_cell{ii} = Simulation.addStar(spectrograph.maxR,spectrograph.pixSamp, simulation.scale, ...
                            star.wavelength,star.spectrum,star.rv);
     
        if curve{ii}.atmosphere == 1
            
            spectral_cell{ii}(:,2) = Simulation.addAtmosphere(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2), ...
                                                           atmosphere.telluric, atmosphere.skyback);
            
        end
        
        spectral_cell{ii}(:,2) = Simulation.addCollectingArea(spectral_cell{ii}(:,2),lbt.apDiameter,lbt.blockFrac);
        spectral_cell{ii}(:,2) = Star.fluxDenToflux(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2));
        
        if curve{ii}.AO == 1 
            
            % if we want ao, first check if its already been made
            
            if exist('AO_throughput','var') == 0
                    
                % if not, make it
                [AO_throughput] = Simulation.combineImagerThroughput(AO_list);
            end
            
            % trim source to WFS band (could change between traces) 
            [wfsWave,wfsSpec] = Simulation.trimToBand(spectral_cell{ii}(:,1),spectral_cell{ii}(:,2),lbti_ao.bandPass*1e-3);
            
            % compute counts on wfs. 
            AOFlux = Simulation.resampleToGrid(AO_throughput(:,1)*1e-3,AO_throughput(:,2),wfsWave);
            wfsCounts = sum(wfsSpec.*AOFlux);
            strehlR = Simulation.calculateStrehlRatio(aoType,seeing,wfsCounts,zenith);
            
        end
        
        if isempty(curve{ii}.throughput) == 0  
            
            if exist('starThroughput','var') ==0 
            
              [starThroughput,tputProg] = Simulation.combineImagerThroughput(star_components);
                
              throughputGrid = Simulation.resampleToGrid(starThroughput(:,1)*1e-3,starThroughput(:,2),spectral_cell{ii}(:,1));
                
              spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2).*throughputGrid;

            end
            
        end
    
        if any(strcmp('fiberLink', curve{ii}.throughput)) == 1
             
            strehl = strehlR; strehlR(:,1) = strehlR(:,1)*1000;% need to work in nm for fiber coupling
            wfe = [0,0,0,0,0,0,0,0]; % complicated input (user can specify wave front error)
            adc = 0; % zenith angle for ads 0-60 in steps of 5
            fiberpos = [0,0,0]; % global position offset in microns (x,y,z)
            dof = 0; % depth of focus (not sure if used yet)
            bandPass = 965:1305; % typical bandpass for fiber coupling
            smfLink = FiberLink(wfe,adc,fiberpos,dof,strehlR,bandPass);
            
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
            
        end
        
        
        % cross disperse spectrum into nOrders orders, trim down wavelengths  
        spectral_cell{ii} = Simulation.Xdisperse(spectral_cell{ii},nOrders,wave_coeff);
        
        if any(strcmp('spectrograph', curve{ii}.throughput)) == 1
            
            % include throughput of spectrograph.
            
            spectral_cell{ii} = Simulation.addSpecThroughput(spectral_cell{ii},spectrograph.finalThroughput,nOrders);
            
            for ii = 1:nOrders
                new_y = Simulation.resampleToGrid(spectrograph.finalThroughput{2}(:,ii),spectrograph.finalThroughput{1}(:,ii),tputProg{1,end}(:,1));
                tempTput{1}(:,ii) = tputProg{1,end}(:,1);
                tempTput{2}(:,ii) = new_y .* tputProg{1,end}(:,2);
            end
            
            %[tempTput] = Simulation.addSpecThroughput(,spectrograph.finalThroughput,nOrders);
            specCell{1} = tempTput;
            specCell{2,1} = spectrograph.name;
            tputProg = [tputProg specCell];
            clear specCell tempTput new_y
            
        end
        
        
%         1) always split into spectral orders - basically copy the spectrum n_order times 
%         3) trim to each spectral order 
%         
%         if we want the throughput then
%             include spectrograph throughput for each order
%         end 
        
% split into spectral orders 
%         
% if throughput is needed  
    % add it
% end        

    elseif strcmp(curve{ii}.source,'etalon')
        %
        %         if throughput
        %
        %             combine them
        %         end
        spectral_cell{ii} = spectral_cell{1};
        
        %
    elseif strcmp(curve{ii}.source,'flat')
        
        
        spectral_cell{ii} = spectral_cell{1};
        %         if throughput
        %             combine them
        %         end
        %     end
    end
    %
end