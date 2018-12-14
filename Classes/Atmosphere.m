classdef Atmosphere
    properties
        telluric    % telluric spectrum (wavelength,flux)
        skyback 	% sky background (wavelength,flux)
        fiberProj   % fiber fiber area on sky
        fiberDiam   % angular fiber diameter 
        scaleFactor % sky background scaling factor
    end
    
    methods
        
        %---------------%
        % Constructor
        %---------------%
        
        function obj = Atmosphere(fiberDiam,scaleFactor)
            
            %creates an atmosphere with telluric absorbtion lines,
            %skybackground scaled to the collecting aperture and assigns
            %an airmass
            
            if nargin == 0
                
                fiberDiam = 43; % default ilocater fiber diameter
                
                scaleFactor = 36; % 36 is scaling factor from gemini to LBT H band.
                
            end
            
            %----------------------------%
            % Load telluric spectrum
            %----------------------------%
            
            current_path = pwd;
            if strcmp(current_path(1:7),'Volumes')==1
            
                telluric_path = 'Volumes/Software/Simulator/RefFiles/Atmosphere/telluric_200.mat';
                skyback_path = 'Volumes/Software/Simulator/RefFiles/SkyBackground/SkyBackground.mat';
            
            elseif strcmp(current_path(2:4),'afs')==1
                telluric_path = '/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/RefFiles/Atmosphere/telluric_200.mat';
                skyback_path = '/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/RefFiles/SkyBackground/SkyBackground.mat';
                
            
            else
                
                telluric_path = [current_path(1:2) '\Simulator\RefFiles\Atmosphere\telluric_200.mat'];
                skyback_path = [current_path(1:2) '\Simulator\RefFiles\SkyBackground\SkyBackground.mat'];

            end
            
            temp = load(telluric_path,'-mat');
            
            temp = struct2cell(temp);
            
            tell = temp{1};
            
            %-----------------------------------%
            % Calculate sky background spectrum
            %-----------------------------------%
            
            temp = load(skyback_path,'-mat');
            
            temp = struct2cell(temp);
            
            skyback = temp;
            
            %-------------------------%
            % Assign object parameters
            %-------------------------%
            
            obj.telluric(:,1)=tell(:,1)/1000; % native telluric file is in nm
            
            obj.telluric(:,2)=tell(:,2); % transmission spectrum
            
            obj.fiberDiam = fiberDiam;
            
            obj.scaleFactor = scaleFactor;
            
            obj.skyback = skyback{1}; % sky background units  = photons/s/nm/m^2/arcsec^2 at LBT magnitude
            
            obj.fiberProj = pi*(0.5*fiberDiam*1E-3)^2; %40mas circle in arcseconds^2
            
            obj.skyback(:,2) = obj.skyback(:,2) * obj.scaleFactor; % convert to LBT from Gemini H band
            
            obj.skyback(:,2) = obj.skyback(:,2)*obj.fiberProj; %Fiber projection on sky will remove /arcsecond^2 in Background
            
            obj.skyback(:,2) = obj.skyback(:,2) * 1000;%nm/mircon converts flux density to microns
            
            obj.skyback(:,1) = obj.skyback(:,1)./1000; % wavelength conversion nm to microns

        end
        
        function [] = plotAtm(obj)
            
            figure
            yyaxis left
            plot(obj.telluric(:,1),obj.telluric(:,2))
            xlabel('wavelength (\mum)')
            ylabel('transmission')
            yyaxis right
            semilogy(obj.skyback(:,1),obj.skyback(:,2),'Color',[0 0 0 0.3])
            hold on
            ylabel('log(counts)')
            box on
            ax=gca;
            ax.LineWidth = 2;
            ax.FontSize = 15;
            ax.YAxis(2).Color = 'k';
            grid on
        end
        
    end
    
end

