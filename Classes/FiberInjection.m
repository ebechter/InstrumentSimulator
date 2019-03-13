classdef FiberInjection < FiberCouplingV2
    
    properties
        name
        Progress
        Fresnel
        Strehl
        bandPass
        StaticSR
        finalThroughput
    end
    
    methods
        function [obj] = FiberInjection(wfe,adc,fiberpos,dof,strehl,bandPass)
            
            %-----------------------%
            % Superclass Constructor
            %-----------------------%
            wavelength =flipud([1.29882 1.28752 1.27643 1.26552 1.25479 1.24425 1.23388 1.22368 1.21365 1.20378 1.19408 ...
                1.18452 1.17512 1.16587 1.15676 1.14779 1.13896 1.13027 1.12171 1.11327 1.10497 1.09678 1.08872 1.08077 ...
                1.07294 1.06522 1.05761 1.05011 1.04271 1.03542 1.02823 1.02114 1.01415 1.00725 1.00044 9.9373E-001 ...
                9.8710E-001 9.8057E-1 9.7411E-1]');
            
            dispersion = true;
            
            %----------Mode Mismatch/ Alignment----------% 
            obj@FiberCouplingV2(wavelength,wfe,adc,fiberpos,dof,dispersion); % call superclass constructor with specific inputs
            
            %----------Assign/ extrapolate parameters----------% 
            obj.name = 'FiberCoupling';
            
            obj.bandPass = bandPass;
            
            obj.Strehl = interp1(strehl(:,1),strehl(:,2),obj.bandPass','linear','extrap');
            
            rho = interp1(obj.Wavelength*1000,obj.Rho,obj.bandPass','linear','extrap'); %mode mismatch overlap integral
            obj.Rho = rho;
            
            %----------Static losses at Fiber----------%         
            obj.Fresnel = 0.92; %fresnel refelction at an uncoated fiber tip and exit      
            
            n1 = 4; n2 = 4;
               
            sig_rms = sqrt(n1*(1/(10*3.5))^2+n2*(1/(20*3.5))^2);
            
            obj.StaticSR = exp(-(2*pi*sig_rms')^2); %Static Strehl ratio delivered by the optics
                       
            %----------Multiply throughput terms----------%
            obj.Progress(:,1) = obj.bandPass';
            obj.Progress(:,2) = obj.Strehl;
            obj.Progress(:,3) = obj.Progress(:,2).*obj.StaticSR;
            obj.Progress(:,4) = obj.Progress(:,3).*obj.Rho;
            obj.Progress(:,5) = obj.Progress(:,4).*obj.Fresnel;
            
            %----------Assign Object values for final throughput----------%
            
            obj.finalThroughput(:,1) = obj.Progress(:,1);
            % obj.finalThroughput(:,2) = ones(size(obj.Progress(:,4),1),1);
            obj.finalThroughput(:,2) = obj.Progress(:,5);
                      
            
        end
        
        
    end
    
end

