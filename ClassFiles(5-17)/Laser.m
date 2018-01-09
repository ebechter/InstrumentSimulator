classdef Laser < Spectra
    properties
        % will fill these out as this class is better understood         
    end
    
    methods
        function [obj] = LoadLaser(obj,laser)
%             if isempty(varagin)==0
%                 laser2 = varagin{1};
%             end
            obj.Wavelength = (laser(1)-0.005: 1e-5 :laser(1)+0.005)' ; % Microns       
            obj.Spectrum = laser(2)*normpdf(obj.Wavelength,laser(1),1e-3);
            
            
            
        end
        
       function [Counts3] = Scale2Star (Counts1, Counts2)
       
           A = max(max(Counts1)); %max of the Etalon spectrum;
           B = max(max(Counts2)); %max of the Stellar spectrum
           ScaleFactor = B/A; %What is the ratio of the stat to cal
           Counts3 = Counts1*ScaleFactor; %Scale the calibration source to the star  
       end
         
    end
end

