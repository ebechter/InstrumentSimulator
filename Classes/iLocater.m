classdef iLocater < Instrument
    properties
    end
    
    methods
        function [obj] = iLocater(entwindow,band,StarlightPath,coupling,wavelength,polarization)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
            else
                % need to put conditions for user specified inputs here
            end
            obj.FullBand = band';
            obj.Polarization = polarization;
            obj = LBTMirrors(obj);
            obj = EntranceWindow(obj,entwindow);
            obj = LBTIOptics(obj);
            obj = CommonChOptics(obj);
            obj = FiberChOptics(obj);
            obj = ImageChOptics(obj);
            obj = WFCOptics(obj);
            obj = QuadOptics(obj);
            obj = SpecOptics(obj);
            obj = R6Grating(obj);
            obj = iLocFiberLink(obj,coupling,wavelength); % calculates fiber coupling
            obj.StrehlRatio = ones(length(obj.FullBand),1);
                        
            %Propogate starlight through the Instrument
            obj.PathName = StarlightPath; %Name of optical path (user input at top)
           
%             obj = Path_Multiply(obj); %Mulitply curves in optical path
            
        end
    end
end