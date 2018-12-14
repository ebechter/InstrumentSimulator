classdef FiberLink
    
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
        bandPass
        specLink
        fwrLink
    end
    
    methods
        function [obj] = FiberLink(bandPass)
            
            %----------Assign/ extrapolate parameters----------% 
            obj.name = 'FiberLink';
            
            obj.bandPass = bandPass';

            %----------Fiber propogation losses----------%
            n = ones(size(bandPass))';
            Avim_connector = 0.98;
            number = 3;
            obj.Avim = Avim_connector^number.*n;
            obj.specSplit = 0.985.*n;
            obj.fwrSplit =  Avim_connector*0.015.*n;
            obj.dBLoss = (1-(10^(1.2/10)*45/1000)).*n; % fibercore spec sheet loss per m
           
            %----------Multiply throughput terms----------%
            obj.Progress(:,1) = obj.bandPass;
            obj.Progress(:,2) = obj.Avim;
            obj.Progress(:,3) = obj.Progress(:,2).*obj.dBLoss;
            obj.Progress(:,4) = obj.Progress(:,3).*obj.specSplit;
            
            %----------Assign Object values for spectrograph link----------%
            
            obj.specLink(:,1) = obj.Progress(:,1);
            obj.specLink(:,2) = obj.Progress(:,4);
            
            %----------Assign Object values for fwr----------%
            obj.fwrLink(:,1) = obj.Progress(:,1);
            obj.fwrLink(:,2) = obj.Progress(:,4).*obj.fwrSplit;
            
            
        end
        
        
    end
    
end

