classdef Atmosphere 
    properties
        telluric
        skyback
        fiberProj
        fiberDiam
        scaleFactor
    end
    methods
        
        %constructor
        function obj = Atmosphere(fiberDiam,scaleFactor)
            if nargin == 0
                fiberDiam = 43; % default ilocater fiber diameter
                scaleFactor = 36; % 36 is scaling factor from gemini to LBT H band.
                
            end
            
                obj.fiberDiam = fiberDiam;
                obj.scaleFactor = scaleFactor; 
                
                pathprefix = pwd;
                telluric_path = '/RefFiles/Atmosphere/telluric_200.mat';
                fullpath = [pathprefix telluric_path];
                temp = load(fullpath,'-mat');
                temp = struct2cell(temp);
                tell = temp{1};
                
                obj.telluric(:,1)=tell(:,1)/1000; % native telluric file is in nm
                obj.telluric(:,2)=tell(:,2); % transmission spectrum
                
                skyback_path = '/RefFiles/SkyBackground/SkyBackground.mat';
                fullpath = [pathprefix skyback_path];
                temp = load(fullpath,'-mat');
                temp = struct2cell(temp);
                obj.skyback = temp{1};
                obj.fiberProj = pi*(0.5*fiberDiam*1E-3)^2; %40mas circle in arcseconds^2
                obj.skyback(:,2) = obj.skyback(:,2) * obj.scaleFactor; % convert to LBT from Gemini H band 
                
                % photons/sec/nm/m^2 at LBT magnitude
                obj.skyback(:,2) = obj.skyback(:,2)*obj.fiberProj; %Fiber projection on sky will remove /arcsecond^2 in Background
                obj.skyback(:,2) = obj.skyback(:,2) * 1000;%nm/mircon
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

