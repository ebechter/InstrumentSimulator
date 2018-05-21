classdef SuperK < Spectra
    properties
        % will fill these out as this class is better understood
    end
    methods
        function [obj] = SuperK()
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
            else
                % need to put conditions for user specified inputs here
            end
            
            obj = LoadSpectrum(obj);
        end
        function [obj] = LoadSpectrum(obj)
            
            path = pwd;
            load(strcat(path,'\RefFiles\Curves\WLS\whiteLase.mat'));
            x(:,1) = whiteLase(:,1);
            y(:,1) = whiteLase(:,2);
            xq = 600:1:1800;
            vq=interp1(x,y,xq);
            
            obj.wavelength(:,1) = xq/1000 ; %nm in native file. convert to in microns
            obj.spectrum(:,1) = 1000*(10.^(vq/10)./1000)  ; %dbm/nm in native file. convert to W/um 
            obj.counts(:,1) = Spectra.energy2Counts(obj.wavelength,obj.spectrum); % convert energy to counts
        end
    end       
end

