classdef Spectra
    properties
        Spectrum
        Counts
        Wavelength
    end
    methods
        function [obj] = Energy2Counts(obj)
            l = obj.Wavelength*1E-6; %m
            c = 2.9979245800 * 10^8; %m/s
            h = 6.62607004E-34; %J*s
            engy = h*c./l;%energy of a photon in Joules
            
            obj.Counts(:,1) = (obj.Spectrum(:,1))./(engy);
        end
        
        function [DsWavelength] = DopplerShift(obj)
            % DESCRIPTION: Doppler shift a spectrum
            c = 2.9979245800 * 10^8;  % Speed of light [m/s] according to NIST - http://physics.nist.gov/cgi-bin/cuu/Value?c
            beta = obj.RV / c;
            delta = sqrt((1 + beta) / (1 - beta));
            DsWavelength = [obj.Wavelength].*delta;
        end
    end
        methods(Static)
            
        function [energy] = Counts2Energy(wavelength,inputspectrum)
            l = wavelength*1E-6; %m
            c = 2.9979245800 * 10^8; %m/s
            h = 6.62607004E-34; %J*s
            engy = h*c./l;%energy of a photon in Joules
            energy = engy.*inputspectrum;
        end
    end
end
