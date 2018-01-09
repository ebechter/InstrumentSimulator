classdef Calibration < Instrument
    properties
    end
    
    methods
        function [obj] = Calibration(band)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                band = (900:1350);
            else
                % need to put conditions for user specified inputs here
            end
            obj.FullBand = band';
            obj = obj.SpecOptics;
            obj = obj.R6Grating;
            
            %Propogate starlight through the Instrument
            obj.PathName = 'Calibration';
            obj = Optical_Path(obj); %Set Optical Path parameters
            obj = Path_Multiply(obj); %Mulitply curves in optical path
            
        end
    end
end