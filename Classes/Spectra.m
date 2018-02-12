classdef Spectra
    properties
        spectrum
        counts
        wavelength
        dsWavelength
        spectrumUnits
    end
    methods
        function [obj] = energy2Counts(obj)
            l = obj.wavelength*1E-6; %m
            c = 2.9979245800 * 10^8; %m/s
            h = 6.62607004E-34; %J*s
            engy = h*c./l;%energy of a photon in Joules
            
            obj.counts(:,1) = (obj.spectrum(:,1))./(engy);
        end
        
    end
    methods(Static)
        
        function [dsWavelength] = dopplerShift(wavelength,rv)
            % DESCRIPTION: Doppler shift a spectrum
            c = 2.9979245800 * 10^8;  % Speed of light [m/s] according to NIST - http://physics.nist.gov/cgi-bin/cuu/Value?c
            beta = rv / c;
            delta = sqrt((1 + beta) / (1 - beta));
            dsWavelength = [wavelength].*delta;
        end
        function [energy] = counts2Energy(wavelength,inputspectrum)
            l = wavelength*1E-6; %m
            c = 2.9979245800 * 10^8; %m/s
            h = 6.62607004E-34; %J*s
            engy = h*c./l;%energy of a photon in Joules
            energy = engy.*inputspectrum;
        end
    end
end
