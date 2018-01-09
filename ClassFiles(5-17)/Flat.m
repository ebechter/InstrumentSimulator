classdef Flat < Spectra
    properties
        % will fill these out as this class is better understood
    end
    
    methods
        
        function [obj] = Flat(scale)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
            else
                % need to put conditions for user specified inputs here
            end
            
            obj = LoadFlat(obj,scale);
            
        end%constructor
        function [obj] = LoadFlat(obj,scale)
            
            R = 275e3;
            pix_samp = 3;
            dlam = 1/R/pix_samp; % average delta lambda per order with 180000 and 3.5 pixel samp.
            step=(dlam/scale);
            
            obj.Wavelength = (0.2: step :1.5)' ; % microns
            obj.Counts = 5e3*ones(size(obj.Wavelength)) ; % BS scale factor
            
        end%generate standard flat
    end
    methods(Static)
        function [Counts3] = Scale2Star (Counts1, Counts2) 
            A = max(max(Counts1)); %max of the Etalon spectrum;
            B = max(max(Counts2)); %max of the Stellar spectrum
            ScaleFactor = B/A; %What is the ratio of the stat to cal
            Counts3 = Counts1*ScaleFactor; %Scale the calibration source to the star
            
        end
        
    end
end

