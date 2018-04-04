classdef Etalon < Spectra
    properties
        FSR        % assigned free spectral range in GHz
        R1         % mirror 1 reflectivity (fractional)
        R2         % mirror 2 reflectivity (fractional)
        finesse    % calculated from R1 and R2.
        length     % cavity length in meters - calculated from FSR
        type       % default, ref, or fullref 
        
        wavecenters % vector of peak wavelengths in angstroms. Calculated by maximizing intensity function  
        refpeak
        %         FSR_measured
        %         FSR_estimated
    end
    
    methods
        function [obj] = Etalon(scale,FSR, R1, R2, rv, type)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                
                FSR = 10;
                R1 = 0.93;
                R2 = 0.93;
                scale = 3;
                rv = 0;
                type = 'default';
            
                
            elseif nargin == 1
                FSR = 10;
                R1 = 0.93;
                R2 = 0.93;
                rv = 0;
                type = 'default';
                
            elseif nargin ==5
                
                type = 'default';
                
            end
            
            obj.FSR = FSR;
            obj.R1 = R1;
            obj.R2 = R2;
            obj.scale = scale;
            obj.rv = rv;
            obj.type = type;
            obj = LoadEtalon(obj); %Andrew added this 5/23/17
            
            if rv ~=0
                obj.rv = rv;
                obj.wavelength = obj.dopplerShift;
            end
        end
        
        
        function [obj] = LoadEtalon(obj)
            
            R1 = obj.R1;
            R2 = obj.R2;
            scale = obj.scale;
            FSR = obj.FSR;
            type = obj.type;
            
            %             [wavelength,Intesnity,p,cent]
            
            %-----User parameters-----%
            
            % lambda = (900:0.0001:1350)*1e-9;
            % scale =6;
            R = 275e3;
            pix_samp = 3;
            
            dlam = 1000/R/pix_samp; % average delta lambda per order with 180000 and 3.5 pixel samp.
            step=(dlam/scale);
            
            lambda = (900:step:1350)*1e-9; % wavelength vector
            
            %-----Constants-----%
            % c = 3.0e8; %Speed of light
            c = physconst('lightspeed');
            n = 1; %refractive index of cavity space
            phi = 0;% 180 degree phase change when reflecting off gold from air
            % FSR = c/(2*l*n)/1E9;%Free spectral range
            obj.length = 1/(FSR*1E9*(2*n/c));
            theta = 0; %normal incidence
            
            %-----Calculations-----%
            delta = (2*pi./lambda)*2*n*obj.length*cos(theta)+2*phi; % ck/FSR
            f_coeff = 4*sqrt(R1*R2)/((1-R1)*(1-R2));
            obj.finesse = pi*(R1*R2)^(1/4)/(1-sqrt(R1*R2));
            % f_coeff_check = (2*finesse/pi).^2;
            I = 1 ./ (1 + f_coeff * sin(delta/2).^2);
            
            iBound = round(2*n*obj.length./[lambda(1) lambda(end)]);
            i = iBound(1):-1:iBound(end);
            obj.wavecenters = 1e10*(2*n*obj.length)./i;
            
            
            
            if strcmp(type, 'default')==1
                cent = 0;
            elseif strcmp(type, 'ref')==1
                %Reference etalon
                lb = 975;
                ub = lb+0.02;
                lambda_ref = (lb:0.000001:ub)*1e-9;
                delta_ref = (2*pi*2*n*obj.length./(lambda_ref)+2*phi);
                I_ref = 1 ./ (1 + f_coeff * sin(delta_ref/2).^2);
                lambda_ref = lambda_ref*1e9; % to nm
                
                [~,ind] = max(I_ref);
                cent = lambda_ref(ind);
                
                Iref1 = I_ref; %single reference etalon peak
                lambda_ref1 = lambda_ref*1e-9; %reference etalon wavelength
                Iref2 = interp1(lambda_ref1,Iref1,lambda,'linear','extrap'); %interpolate reference onto science band
                Iref3 = Iref2*5; %scale up reference peak
                Ind = (lambda>lb*1e-9 &lambda<ub*1e-9); %original band = 1, outside band = 0
                Iref4 = Iref3.*Ind; %mutliplty logical by Intensity to set outside band to 0 (extrap doesnt work well)
                Iref5 = Iref4+I; %add bright peak to science etalon
                I = Iref5;% overwrite the original with a bright peak
                
                
            elseif strcmp(type, 'fullref')==1
                %Reference etalon
                
                startpoints = [0.974776348372602,0.981232163872611,0.987774054151196,0.994403757468148,1.00112305300997,1.00793376679567,1.01483778074376,1.02183702204521,1.02893347272005,1.03612917785920,1.04342623213786,1.05082679095501,1.05833307073641,1.06594735823831,1.07367199173734,1.08150939451799,1.08946206094173,1.09753253890670,1.10572347541185,1.11403758257728,1.12247764059590,1.13104657470581,1.13974731358388,1.14858296767413,1.15755670296739,1.16667173064714,1.17593144228870,1.18533931447052,1.19489893770972,1.20461399498013,1.21448828040811,1.22452577664988,1.23473060388879,1.24510688571311,1.25565906336888,1.26639164636785,1.27730920927658,1.28841671281211,1.29971899892678];
                [pks,locs] = findpeaks(I,lambda);
                
                
                for ii = 1:39
                    [~,ind_peak]= min(abs(locs-startpoints(ii)*1e-6));
                    %         [bs]= find(locs > startpoints(ii)*1e-6,1,'first');
                    test(ii) = ind_peak-bs;
                    roughpeak = locs(ind_peak);
                    lb = roughpeak*1e9 - 0.01;
                    ub = roughpeak*1e9 + 0.01;
                    
                    lambda_ref = (lb:0.000001:ub)*1e-9;
                    delta_ref = (2*pi*2*n*length./(lambda_ref)+2*phi);
                    I_ref = 1 ./ (1 + f_coeff * sin(delta_ref/2).^2);
                    lambda_ref = lambda_ref*1e9; % to nm
                    
                    [~,ind] = max(I_ref);
                    cent(ii) = lambda_ref(ind);
                    
                    Iref1 = I_ref; %single reference etalon peak
                    lambda_ref1 = lambda_ref*1e-9; %reference etalon wavelength
                    Iref2 = interp1(lambda_ref1,Iref1,lambda,'linear','extrap'); %interpolate reference onto science band
                    Iref3 = Iref2*5; %scale up reference peak
                    Ind = (lambda>lb*1e-9 &lambda<ub*1e-9); %original band = 1, outside band = 0
                    Iref4 = Iref3.*Ind; %mutliplty logical by Intensity to set outside band to 0 (extrap doesnt work well)
                    Iref5(ii,:) = Iref4;
                end
                Iref6 = sum(Iref5,1);
                Iref7 = Iref6+I; %add bright peak to science etalon
                I = Iref7;% overwrite the original with a bright peak
            end
            
            
            %-----Assign values to Etalon object-----%
            obj.wavelength(:,1) = lambda*1e6; % microns
            obj.counts(:,1) = I; % BS scale factor
            
            obj.refpeak = cent;
        end
        
        function [] = PlotEtalon (obj)
            close all
            %             figure()
            %             plot(obj.Wavelength/10,obj.dlambda)
            
            figure()
            plot(obj.wavelength/10, obj.counts , 'linewidth' , 0.5 , 'color' , [0.8 , 0.1 ,0.8]);
            ylim([0 1.2])
            
            % %             window = 5*0.1*1E9;
            %             center = mean(obj.Wavelength)/10;
            %             lb = center-window;
            %             ub = center+window;
            % %             xlim([lb ub])
            
            disp(['FSR = ',num2str(obj.FSR),'GhZ'])
            disp(['F = ',num2str(obj.Finesse)])
        end
        
        
        % Potato calculation of FSR.
        function [obj] = FSRCalc(obj)
            
            I1 = obj.Counts;
            lambda = obj.Wavelength*1e-10;
            d = obj.Length;
            
            [pks, locs] = findpeaks(I1);
            nsamps = 20;
            
            a1 = 0.5e-3;
            b1 = lambda(locs(1)); %reset every loop iteration
            c1 = 0;
            d1 = 0.02;
            
            lorentzian = 'd1 /(pi)*((0.5*a1)./((x-b1).^2 + (0.5*a1^2)))+c1';
            startpoints = [a1 b1 c1 d1];
            
            tic
            
            for ii = 1:2;%length(pks)
                
                startpoints(2) = lambda(locs(ii))*1e9; %into nm
                ydata  = I1(locs(ii)-nsamps:locs(ii)+nsamps);
                xdata = lambda(locs(ii)-nsamps:locs(ii)+nsamps)*1e9;
                %     figure(100)
                %     plot(xdata,ydata,'ok')
                myfit2 = fit(xdata,ydata,lorentzian,'Start',startpoints);
                params2 = coeffvalues(myfit2);
                centers(ii) = params2(2);
                %     hold on
                %     plot(myfit2)
            end
            toc
            
            obj.FSR_measured = mean(diff(centers));
            
            lambda0 = centers(1)*1e-9;
            peakcenter = lambda0;
            ii = 1;
            while lambda0 <= lambda(end)% compare in m
                
                dlambda = (lambda0.^2)./(2*1*d); % calc in m
                lambda0 = lambda0+dlambda;
                peakcenter = [peakcenter lambda0];
                FSR_temp(ii) = dlambda;
                ii = ii+1;
            end
            
            obj.FSR_estimated =(diff(peakcenter));
            
        end
        
    end
    methods(Static)
        
    end
end







