classdef FiberLink < FiberCoupling
    
    properties
        name
        fiberLength
        numInterfaces
        splittingRatio
        fresnelLoss
        Progress
        specSplit
        fwrSplit
        Fresnel
        dBLoss
        Avim
        Strehl
        bandPass
        specLink
        fwrLink
    end
    
    methods
        function [obj] = FiberLink(wfe,adc,fiberpos,dof,strehl,bandPass)
         obj.name = 'fiberLink';
         % similar constructor to the fibercoupling superclass but with
         % specific ilocater defaults and static losses through the
         % entire fiberLink

            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                obj.WFE = [0,0,0,0,0,0,0,0]; % Pupil Zernike amplitudes piston, xtilt, ytilt, etc...
                rms_scale = [1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8),1/sqrt(8)];
                obj.RMS = obj.WFE.*rms_scale;
                obj.ADC = 0; %Choose values from 0-60 in steps of 5 for the zenith angle.(This is translated as an offset fiber pos rather than beam position for computational ease)
                obj.FiberPos = [0,0,0]; % global position offset in microns (x,y,z)
                obj.DoF = 0; %Polychromatic depth of focus in microns
                obj.bandPass = 965:1305;
%                 obj.Strehl(:,1) = bandPass;
                obj.Strehl = ones(size(obj.bandPass,2),1);
                
            else
                obj.WFE = wfe; % Pupil Zernike amplitudes piston, xtilt, ytilt, etc...
                rms_scale = [1,1/2,1/2,1/sqrt(3),1/sqrt(6),1/sqrt(6),1/sqrt(8),1/sqrt(8)];
                obj.RMS = obj.WFE.*rms_scale;                
                obj.ADC = adc; %Choose values from 0-60 in steps of 5 for the zenith angle.(This is translated as an offset fiber pos rather than beam position for computational ease)
                obj.FiberPos = fiberpos; % global position offset in microns (x,y,z)
                obj.DoF = dof; %Polychromatic depth of focus in microns
                obj.bandPass = bandPass;
                obj.Strehl = interp1(strehl(:,1),strehl(:,2),obj.bandPass','linear','extrap');

            end
            
            obj.Wavelength =flipud([1.29882 1.28752 1.27643 1.26552 1.25479 1.24425 1.23388 1.22368 1.21365 1.20378 1.19408 ...
                1.18452 1.17512 1.16587 1.15676 1.14779 1.13896 1.13027 1.12171 1.11327 1.10497 1.09678 1.08872 1.08077 ...
                1.07294 1.06522 1.05761 1.05011 1.04271 1.03542 1.02823 1.02114 1.01415 1.00725 1.00044 9.9373E-001 ...
                9.8710E-001 9.8057E-1 9.7411E-1]');
            
            %---------chromatic Losses----------%
            %run fiber coupling with default inputs
            
            obj = PrepareFiberOffsets(obj);
            obj = CreatePupilPlane(obj);
            obj = CreateFiberPSF(obj);
            obj = CoupleSMF(obj);
            
            rho = interp1(obj.Wavelength*1000,obj.Rho,obj.bandPass','linear','extrap'); %mode mismatch overlap integral
            obj.Rho = rho;
            
            %----------Light Losses realted to Fiber injection parameters----------%         
            obj.Fresnel = 0.92; %fresnel refelction at an uncoated fiber tip and exit      
            
            %Static Strehl ratio delivered by the optics
%             [SR_Static]=Wagner();
            
            %----------Light Losses realted to Fiber propogation----------%
            Avim_connector = 0.98;
            number = 3;
            obj.Avim = Avim_connector^number;
            obj.specSplit = obj.Avim*0.99;
            obj.fwrSplit = obj.Avim*0.01;
            obj.dBLoss = 1-(10^(1.2/10)*45/1000); % fibercore spec sheet loss per m
           
            %----------Multiply throughput terms----------%
            obj.Progress(:,1) = obj.bandPass';
            obj.Progress(:,2) = obj.Strehl;
            obj.Progress(:,3) = obj.Progress(:,2).*obj.Rho;
            obj.Progress(:,4) = obj.Progress(:,3).*obj.Fresnel;
            obj.Progress(:,5) = obj.Progress(:,4).*obj.Avim;
            obj.Progress(:,6) = obj.Progress(:,5).*obj.dBLoss;
            obj.Progress(:,7) = obj.Progress(:,6).*obj.specSplit;
           
            %----------Assign Object values for spectrograph link----------%
            
            obj.specLink(:,1) = obj.Progress(:,1);
            obj.specLink(:,2) = obj.Progress(:,7);
            
            %----------Assign Object values for fwr----------%
            obj.fwrLink(:,1) = obj.Progress(:,1);
            obj.fwrLink(:,2) = obj.Progress(:,6).*obj.fwrSplit;
            
            
        end
        
        
    end
    
end

