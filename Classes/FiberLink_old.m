classdef FiberLink_old < FiberCoupling
    
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
        StaticSR
    end
    
    methods
        function [obj] = FiberLink_old(wfe,adc,fiberpos,dof,strehl,bandPass)
            
            %-----------------------%
            % Superclass Constructor
            %-----------------------%
            
            %----------Mode Mismatch/ Alignment----------% 
            obj@FiberCoupling(wfe,adc,fiberpos,dof); % call superclass constructor with specific inputs
            
            %----------Assign/ extrapolate parameters----------% 
            obj.name = 'FiberLink_old';
            
            obj.bandPass = bandPass;
            
            obj.Strehl = interp1(strehl(:,1),strehl(:,2),obj.bandPass','linear','extrap');
            
            rho = interp1(obj.Wavelength*1000,obj.Rho,obj.bandPass','linear','extrap'); %mode mismatch overlap integral
            obj.Rho = rho;
            
            %----------Static losses at Fiber----------%         
            obj.Fresnel = 0.92; %fresnel refelction at an uncoated fiber tip and exit      
            
            n1 = 4; n2 = 4;
               
            sig_rms = sqrt(n1*(1/(10*3.5))^2+n2*(1/(20*3.5))^2);
            
            obj.StaticSR = exp(-(2*pi*sig_rms')^2); %Static Strehl ratio delivered by the optics
            
%             %----------Fiber propogation losses----------%
            Avim_connector = 0.98;
            number = 3;
            obj.Avim = Avim_connector^number;
            obj.specSplit = 0.985;
            obj.fwrSplit =  Avim_connector*0.015;
            obj.dBLoss = 1-(10^(1.2/10)*45/1000); % fibercore spec sheet loss per m
           
            %----------Multiply throughput terms----------%
            obj.Progress(:,1) = obj.bandPass';
            obj.Progress(:,2) = obj.Strehl;
            obj.Progress(:,3) = obj.Progress(:,2).*obj.StaticSR;
            obj.Progress(:,4) = obj.Progress(:,3).*obj.Rho;
            obj.Progress(:,5) = obj.Progress(:,4).*obj.Fresnel;
            obj.Progress(:,6) = obj.Progress(:,5).*obj.Avim;
            obj.Progress(:,7) = obj.Progress(:,6).*obj.dBLoss;
            obj.Progress(:,8) = obj.Progress(:,7).*obj.specSplit;
            
            %----------Assign Object values for spectrograph link----------%
            
            obj.specLink(:,1) = obj.Progress(:,1);
            % obj.specLink(:,2) = ones(size(obj.Progress(:,7),1),1);
            obj.specLink(:,2) = obj.Progress(:,7);
            
            %----------Assign Object values for fwr----------%
            obj.fwrLink(:,1) = obj.Progress(:,1);
            obj.fwrLink(:,2) = obj.Progress(:,6).*obj.fwrSplit;
            
            
        end
        
        
    end
    
end

