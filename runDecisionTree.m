for ii = tracenum
   

    if strcmp(curve{ii}.source,'star')
    spectral_cell{ii}(:,1) = star.wavelength;
    
        if curve{ii}.atmosphere == 1 
            
            spectral_cell{ii}(:,2) = Simulation.addAtmosphere(star.counts,star.wavelength, ...
                                                           atmosphere.telluric, atmosphere.skyback);
            
        end
    
        if curve{ii}.AO == 1 
            
            % if we want ao, first check if its already been made
            
            if exist('AO_throughput','var') ==0
                    
                % if not, make it
                % [AO_throughput] = combinedImagerThroughput(objects);
            else
                add it
                
                % spectral_cell{ii} = Simulation.addAO(spectral_cell{ii},AO_throughput);
            end
        
        end                                                                                                                                                                                                                                                                                                                                                                       
        if isempty(curve{ii}.throughput) == 0  
            
            if exist('starThroughput','var') ==0 
            
              starThroughput = Simulation.combineImagerThroughput(star_components);
                
              throughputGrid = Simulation.resampleToGrid(starThroughput(:,1)*1e-3,starThroughput(:,2),star.dsWavelength);
                
              spectral_cell{ii}(:,2) = spectral_cell{ii}(:,2).*throughputGrid;

            end
            
        end
    
        
        
        % cross disperse spectrum into nOrders orders, trim down wavelengths  
        spectral_cell{ii} = Simulation.Xdisperse(spectral_cell{ii},nOrders,wave_coeff);
        
        if any(strcmp('spectrograph', curve{ii}.throughput)) == 1
            
            % combine throughput of spectrograph. 
            
            spectral_cell{ii} = Simulation.addSpecThroughput(spectral_cell{ii},spectrograph.finalThroughput,nOrders);
            
            
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