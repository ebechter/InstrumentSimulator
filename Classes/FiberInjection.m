classdef FiberInjection < FiberCoupling
    
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
            
            %----------Mode Mismatch/ Alignment----------% 
            obj@FiberCoupling(wfe,adc,fiberpos,dof); % call superclass constructor with specific inputs
            
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

