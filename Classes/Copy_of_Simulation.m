classdef Simulation < Star
    properties
        StrehlRatio
        InstrumentSpectrum
        AOSpectrum
        IntegratedCounts
        IntegratedEnergy
        SpectrumOrders
        EtalonOrders
        FlatOrders
        Airmass
        RawSpectrum
        Bandpass
        DetectorFace
        FWRSpectrum
        SpectrumUnits
        
    end
    
    methods
        %constructor
        function obj = Simulation(type,mag,epsilon,vsini,rv,z_deg,units,fiber_mas)
            if nargin == 0
                obj.SpType = 'M0V'; %spectral type
                obj.Vmag = 10; % V magntiude
                obj.Epsilon = 1; % not sure about this parameter
                obj.Vsini= 2.5; %rotational velocity in km/s (Rot broad uses km/s)
                obj.RV = 0; %Set RV shift
                z_deg = 45; % Zenith angle used to set airmass
                obj.SpectrumUnits = 'counts'; 
                disp('Input parameters may not be set correctly...using defaults')
                obj.SkyBackground = obj.Skyback(40);
            else
                obj.SpType = type; %spectral type
                obj.Vmag = mag; % V magntiude 
                obj.Epsilon = epsilon; % not sure about this parameter
                obj.Vsini= vsini; %rotational velocity in km/s (Rot broad uses km/s)
                obj.RV = rv; %Set RV shift
                obj.SpectrumUnits = units;
                obj = obj.SkyBack(fiber_mas);

            end
            %Exectute default methods on object for all instances. 
            
            %Load the calibrated spectrum. %Scales to magnitude and calculates color to temp
            obj = obj.LoadSpectrum; %target.Spectrum is in w/m^2/micron, target.Wavelength in microns
            %-----Stellar Effects-----%
            obj.Spectrum = obj.RotBroad; %Broaden spectrum and overwrite target.Spectrum property
            %Convert Spectum to counts and fills target.Counts property
            obj = Energy2Counts(obj); %target.Counts is in counts/s/m^2/micron
            obj.DsWavelength = obj.DopplerShift; %Shift Wavelength and assign DsWavelength
            %-----Atmospheric Effects-----%
            obj = obj.LoadTelluric;%Trims DsWavelength, Wavelength, Tellurics, and Spectrum
            
            %Raw spectrum values are set for reference.
            obj.RawSpectrum(:,1) = obj.Wavelength;
            obj.RawSpectrum(:,2) = obj.Spectrum;
            obj.RawSpectrum(:,3) = obj.Counts;
            obj.RawSpectrum(:,4) = obj.Tellurics;
            
            %-----Observation Conditions -----%
            % target = target.SkyBack; %Adds property Skybackground to object trimmed to DsWavelength property
            obj.Airmass(1,1) = z_deg;
            obj.Airmass(1,2) = z_deg*pi/180;
            obj.Airmass(1,3) = sec(obj.Airmass(1,2));
            
            
            
        end
        function obj = CalculateStrehl(obj,FLAO)
            global strehl_dir
            counts=obj.IntegratedCounts;
            
            if FLAO == 1
                filename = 'FLAO_Strehl.csv';
            elseif FLAO==2
                filename = 'SOUL_Strehl_1_vib.csv';
            elseif FLAO==3
                filename = 'SOUL_Strehl_08_vib.csv';
            elseif FLAO==4
                filename = 'SOUL_Strehl_06_vib.csv';
            end
            
            [MagR,SRY,SRJ,SRR,SRI,SRH,SRK,SRL,SRM,fsHz,fluxphfrpix,bin,fluxphfr] = importStrehl(strcat(strehl_dir,filename));
            clear path filename
            
            phfr= (fluxphfr(:,1).*(bin.^2));
            phs = phfr(:,1).*fsHz;
            phfrpix_binned = fluxphfrpix(:,1).*(bin.^2);
            phfrpix = fluxphfrpix(:,1);
            phfrpixppl = phfrpix/4;
            
            high = max(phs);
            low  = 0;
            num = 1000000;
            xq = linspace(low,high,num);
            vqy = interp1(phs,SRY,xq,'pchip');
            vqj = interp1(phs,SRJ,xq,'pchip');
            vqr = interp1(phs,SRR,xq,'pchip');
            vqi = interp1(phs,SRI,xq,'pchip');
            vqh = interp1(phs,SRH,xq,'pchip');
            vqk = interp1(phs,SRK,xq,'pchip');
            vql = interp1(phs,SRL,xq,'pchip');
            vqm = interp1(phs,SRM,xq,'pchip');
            
            if counts > high
                StrehlY = max(vqy);
                StrehlJ = max(vqj);
                StrehlR = max(vqr);
                StrehlI = max(vqi);
                StrehlH = max(vqh);
                StrehlK = max(vqk);
                StrehlL = max(vql);
                StrehlM = max(vqm);
                ind = 1;
            else
                ind = find(xq>counts,1,'first');
                StrehlY = vqy(ind);
                StrehlJ = vqj(ind);
                StrehlR = vqr(ind);
                StrehlI = vqi(ind);
                StrehlH = vqh(ind);
                StrehlK = vqk(ind);
                StrehlL = vql(ind);
                StrehlM = vqm(ind);
            end
            
            if FLAO ==1
                fs = [1000.1 1000.01 1000.001 1000.0001 200 100];
                bin = [1.0001 1.001 2 4 4.0001 4.001 4.01];
            else
                fs = [1500.1 1500.01 1250 750 300 100];
                bin = [1.0001 1.001 1.01 2 3 4];
            end
            
            %Some WFS parameters....
%             upsamp_fsHz = interp1(phs,fs,xq);
%             
%             pfr_counts = xq(ind)/round(upsamp_fsHz(ind));
%             FR = round(upsamp_fsHz(ind));
%             
%             upsamp_bins = interp1(phs,bin,xq);
%             BIN = round(upsamp_bins(ind));
            
            % figure()
            % hold on
            % plot(xq,vqy,'.')
            % plot(phs,SRY,'.')
            % plot(xq,vqj,'.')
            % plot(phs,SRJ,'.')
            % plot(xq,vqr,'.')
            % plot(phs,SRR,'.')
            
            StrehlSamples = [StrehlR,StrehlI,StrehlY,StrehlJ,StrehlH,StrehlK,StrehlL,StrehlM];
            
            [SR_interp]=SR_interpolation(StrehlSamples);
            
            obj.StrehlRatio = SR_interp;
            
        end
        function obj = PlotSpectrumOrders(obj)
            figure
            hold on
            for ii = 1:size(obj.SpectrumOrders{1},2)
                plot(obj.SpectrumOrders{1}(:,ii),obj.SpectrumOrders{2}(:,ii))
            end
            xlim(obj.Bandpass(1,:))
            xlabel('\lambda \mu m')
            box on
        end
        function obj = PlotEtalonOrders(obj)
            figure
            hold on
            for ii = 1:size(obj.EtalonOrders{1},2)
                plot(obj.EtalonOrders{1}(:,ii),obj.EtalonOrders{2}(:,ii))
            end
            xlim(obj.Bandpass(1,:))
            xlabel('\lambda \mu m')
            box on
        end
        function obj = PlotFlatOrders(obj)
            figure
            hold on
            for ii = 1:size(obj.FlatOrders{1},2)
                plot(obj.FlatOrders{1}(:,ii),obj.FlatOrders{2}(:,ii))
            end
            xlim(obj.Bandpass(1,:))
            xlabel('\lambda \mu m')
            box on
        end     
        function obj = PlotInstrumentSpectrum(obj)
            global colors
            figure()
            hold on
            plot(obj.InstrumentSpectrum(:,1),obj.InstrumentSpectrum(:,2),'color',colors{8},'LineWidth',0.01)
            xlabel('\lambda (um)')
            ylabel('Spectral Intensity (Ph s^-1 um^-1)')
            [bandpass,flux_cut]=select_bandpass(obj.InstrumentSpectrum(:,1),obj.InstrumentSpectrum(:,2),obj.Bandpass(2,:));
            plot(bandpass,flux_cut,'color',colors{1});
            yl = ylim;
            ylim([0, yl(1,2)])
            box on
            
        end
        function obj = PlotSpectrum(obj)
            
        end
        function obj = NumericalThroughputGrating(obj)
            
            for ii = 1:size(obj.SpectrumOrders{2},2)
                obj.SpectrumOrders{1,5}(:,ii) = obj.SpectrumOrders{1,2}(:,ii)./obj.SpectrumOrders{1,4}(:,ii);
            end
            
        end
        
    end
    
    methods(Static)
        function [spectrum_new,wavelength_new] = IncludeThroughput...
                (throughput_wavelength,spectrum_wavelength,Tput,spectrum)
            lb = max(throughput_wavelength(1),spectrum_wavelength(1));
            ub = min(throughput_wavelength(end),spectrum_wavelength(end));
            ind2 = find(spectrum_wavelength<=ub,1,'last');
            ind = find(spectrum_wavelength>=lb,1,'first');
            
            spectrum_wavelength_new = spectrum_wavelength(ind:ind2);
            spectrum_new = spectrum(ind:ind2);
            tput_new = interp1(throughput_wavelength(:,1),Tput(:,1),spectrum_wavelength_new(:,1));
%             tput_new = 0.01*ones(size(spectrum_new)); %removes effects of
%             throughput - This should be obsolete now 6/5/17
            wavelength_new = spectrum_wavelength_new;
            spectrum_new =spectrum_new.*tput_new;
            %now use f and wave(um) for flux and wavelength range within the instrument
           
        end
        function [spectrum] = LBT_CollectingArea(spectrum)
            tel_rad = 4.11;
            spectrum = spectrum*(pi*(tel_rad^2));
        end
        function [SR_out] = rescale_SR(SR_in,z)
            %--------------------------------------------------------------------------                                        General
            %Description: Rescale the Strehl ratio as a function of zenith angle
            %Input: Strehl ratio from 0-1, zenith angle in radians
            %Documentation:
            %Roddier A0 Textbook and Devaney 2007
            %SR = exp(-sigma^2)
            %r0 ~ cos(z)^3/5 scale factor for seeing angle (from Francois Roddier AO book)
            %The intention is to rescale the freid parameter which is buried in a
            %constant exp -(A/ro)^-5/3 or exp(-B). ro scales as cos(z)^3/5 then the constant (B) scales
            %as (cos(z)^3/5)^-5/3 which is just sec(z)
            %--------------------------------------------------------------------------
            
            % if max(SR_in)>1 %the following calculation only works in SR 0-1
            SR = SR_in./100;
            % end
            B = abs(log(SR));%determine the constant
            scale = sec(z);%.^(3/5).^(5/3); % also known as sec(z)
            SR_z = exp(-B.*scale);% rescale the strehl ratio and multiply by 100 to get back to % Strehl
            
            % if max(SR_in)>1 % if the orginal strehl was in %, adjust the modified strehl to be in %
            SR_out = SR_z.*100;
            % else
            %     SR_out = SR_z; % in the original strehl was from 0-1, do nothing
            % end
            
        end
        %Read in a Zemax Kernel
        function[F,kernel_prop] = ReadZemaxKernel(filename,varargin)
            
            %% Reads Zemax Kernel
            %% Fixes non-square physical grid dimensions output from Zemax
            %% Rescales to specific sampling, if specified.
            
            %% need to include horizonal xlength value
            
            formatSpec{1} = '%s %f %s %s %s %f %s %f %s %f %s\n';
            formatSpec{2} = '%s %s %s %s %f %s %s %s %s %f %s\n';
            
            fileID = fopen(filename,'r','n','UTF-8');
            A = textscan(fileID,formatSpec{1},'HeaderLines',10);
            B = textscan(fileID,formatSpec{2},'HeaderLines',0);
            fclose(fileID);
            
            % kernel properties
            kernel_prop = [A{2} B{5} B{10} A{8}];
            % [lambda,X,Y] microns, mm,mm, xcen(mm)
            
            % sprintf('Kernel Properties:\n\nlambda = %.5f microns\nX = %.4f mm\nY = %.4f mm\n\n ',kernel_prop)
            
            kernel_array=dlmread(filename,'\t',14,0);
            
%             Check initial array
%             figure
%             imagesc(kernel_array)
            
            
            % Create x and y vectors in mm
            xx=0:kernel_prop(2)/(size(kernel_array,1)-1):kernel_prop(2);
            yy=0:kernel_prop(3)/(size(kernel_array,2)-1):kernel_prop(3);
            
            % always scale to the minimum of x and y
            scaleto = min(kernel_prop(2),kernel_prop(3));
            [X,Y] = meshgrid(xx,yy);
            
            if scaleto==xx(end)
                [X2,Y2] = meshgrid(xx,xx);
                kernel_prop(3)=kernel_prop(2);
            else
                [X2,Y2] = meshgrid(yy,yy);
                kernel_prop(2)=kernel_prop(3);
            end
            
            
            % interpolate to new physically "square" grid
            kernel_array = interp2(X,Y,kernel_array,X2,Y2);
            
            % sprintf('New Kernel Properties:\n\nlambda = %.5f microns\nX = %.4f mm\nY = %.4f mm\n\n ',kernel_prop)
            %% Check if a rescaling parameter is specified
            if isempty(varargin) ~=1
                
                sum_array = sum(kernel_array);
                
                % fit gaussian to x direction
                [F,GOF] = fit((1:size(kernel_array,2))',sum_array','a*exp(-((x-b)./(2.*c)).^2)','Start'...
                    ,[max(max(kernel_array)) size(kernel_array,2)/2 size(kernel_array,2)/5]);
                
                coeffs = coeffvalues(F);
                sigma = 2*sqrt(2*log(2))*coeffs(3);
                
                % determine scaling for new grid
                scale_factor = varargin{1}/sigma;
                
                % resize kernel using scale factor
                a= imresize(kernel_array,scale_factor);
                % add a bias to fix negatives.
                a = a+abs(min(min(a)));
                % normalize area to 1
                a = a./sum(sum(a));
                
                % Testing materials
                [MatX,MatY]=meshgrid(1:1:size(a,1),1:1:size(a,2));     
                F=circ_gauss(MatX,MatY,varargin{1}/2*sqrt(2*log(2)),[mean(1:size(a,1)),mean(1:size(a,1))]);

                
                
                % check
                %     [G,GOF] = fit((1:size(a,2))',sum(a)','a*exp(-((x-b)./(2.*c)).^2)','Start'...
                %         ,[max(max(a)) size(a,2)/2 varargin{1}/2]);
                
            end
            
            
        end
        %Convert wavelength solution information (fits zemax smaples)
        function [tel_data,xy_coeff,lamx_coeff] = ConvertZemaxFormat(file)
            %% Read in spectral format file.
            order_data=dlmread(file);
            order_data=flipud(order_data);
            % 1 Telescope data at a time
            tel1 = read_orders(1,order_data);
            tel2 = read_orders(2,order_data);
            tel3 = read_orders(3,order_data);
            %% XY
            % figure()
            % hold on
            for ii = 1:36
                %     hold on
                %
                %     % fit each order with a polynomial, save the coefficients
                xy_coeff(ii,1:3,1) = polyfit(tel1(ii,:,1),tel1(ii,:,2),2);
                xy_coeff(ii,1:3,2) = polyfit(tel2(ii,:,1),tel2(ii,:,2),2);
                xy_coeff(ii,1:3,3) = polyfit(tel3(ii,:,1),tel3(ii,:,2),2);
                
                % %       Plot all the orders
                %     plot(tel1(ii,:,1),polyval(xy_coeff(ii,1:3,1),tel1(ii,:,1)))
                %     plot(tel2(ii,:,1),polyval(xy_coeff(ii,1:3,2),tel2(ii,:,1)))
                %     plot(tel3(ii,:,1),polyval(xy_coeff(ii,1:3,3),tel3(ii,:,1)))
                %     pause
            end
            %% lam-X
            % figure()
            for ii = 1:36
                %     hold on
                %     plot(tel1(ii,:,1),tel1(ii,:,3),'o')
                
                % fit each order with a polynomial, save the coefficients
                % switching to 4th order 9/6/17
                lamx_coeff(ii,:,1) = polyfit(tel1(ii,:,1),tel1(ii,:,3),4);
                lamx_coeff(ii,:,2) = polyfit(tel2(ii,:,1),tel2(ii,:,3),4);
                lamx_coeff(ii,:,3) = polyfit(tel3(ii,:,1),tel3(ii,:,3),4);
                %
                %       Plot all the orders
                %     plot(tel1(ii,:,1),polyval(lamx_coeff(ii,:,1),tel1(ii,:,1)))
                %     plot(tel2(ii,:,1),polyval(lamx_coeff(ii,:,2),tel2(ii,:,1)))
                %     plot(tel3(ii,:,1),polyval(lamx_coeff(ii,:,3),tel3(ii,:,1)))
            end
            %% X-lam -- using 3rd order polynomial
            % figure()
%             for ii = 1:39
%                 %     hold on
%                 %     plot(tel1(ii,:,3),tel1(ii,:,1),'o')
%                 
%                 % fit each order with a polynomial, save the coefficients
%                 xlam_coeff(ii,:,1) = polyfit(tel1(ii,:,3),tel1(ii,:,1),3);
%                 xlam_coeff(ii,:,2) = polyfit(tel2(ii,:,3),tel2(ii,:,1),3);
%                 xlam_coeff(ii,:,3) = polyfit(tel3(ii,:,3),tel3(ii,:,1),3);
%                 %
%                 % %       Plot all the orders
%                 %     plot(tel1(ii,:,3),polyval(xlam_coeff(ii,:,1),tel1(ii,:,3)))
%                 %     plot(tel2(ii,:,3),polyval(xlam_coeff(ii,1:3,2),tel2(ii,:,3)))
%                 %     plot(tel3(ii,:,3),polyval(xlam_coeff(ii,1:3,3),tel3(ii,:,3)))
%             end
            
            
            
            tel_data(1,:,:,:) = tel1;
            tel_data(2,:,:,:) = tel2;
            tel_data(3,:,:,:) = tel3;
            
        end
        % Convolve orders with Zemax kernel
        function [PSF,center] = MakePSF(scale,pix_samp,vert_samp)
          
            % dispersion direction
            fwhm = pix_samp*scale;
            sigmax=fwhm/(2*sqrt(2*log(2)));
            
            % cross dispersion direction
            fwhmy = vert_samp*scale;
            sigmay=fwhmy/(2*sqrt(2*log(2)));
           
            % trim according to bigger dimension
            Nsigmas = 9;
            val = Nsigmas * sigmay;
%                         
            
%             clip = 3/(2*sqrt(2*log(2)));
%             val = 25*clip;
            
            oddGrid=2*round(val/2)-1; % not super important but lets go with an odd grid size every time
            
            [MatX,MatY]=meshgrid(1:1:oddGrid,1:1:oddGrid); % make the mesh grid
            
            center = [round(size(MatX,2)/2),round(size(MatY,1)/2)]; % center x and y in the (likely) non-integer center of grid ;
            % this is the most important line - the center of the PSF must fall in a pixel center
            % subtract the center in the x direction from the full convolution to get a "central" conv. result
            %....maybe not certain.
            
            PSF=circ_gauss(MatX,MatY,[sigmax,sigmay],center);

            
        end
        function [trim, wavelength] = ConvolveOrder(wavelength,spectrum,wave_coeff,scale)
            
            % trim each order beyond the edge of the detector
            
            ind1 = find(wavelength <= polyval(wave_coeff,-22.48),1,'last');
            ind2 = find(wavelength >= polyval(wave_coeff, 22.48),1,'first');
            
            wavelength = wavelength(ind1:ind2)';
            spectrum = spectrum(ind1:ind2)';
            
            
            % sampled at the high end. 3 pixels at red, smooth function to
            % 6 pixels at blue. 
            
               
            pix_samp = linspace(6,3,length(wavelength));

            vert_samp = median(pix_samp);
            pix_samp = 0.8*pix_samp;            
            % Custom convolution
            
            % Do the first loop iteration outside the loop. Need to
            % calculate dim first
            ii = 1;
            
            [PSF,center] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);           
            dim = size(PSF,1);         
            rectangle = zeros(dim,length(wavelength)+dim-1);            
            rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
            
            for ii = 2:length(wavelength)
                    
                [PSF,~] = Simulation.MakePSF(scale,pix_samp(ii),vert_samp);
                
%                 Conv1 = conv2(spectrum(ii),PSF,'full');
                
                rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+PSF.*spectrum(ii);
%                 rectangle(:,ii:dim-1+ii)=rectangle(:,ii:dim-1+ii)+Conv1;
                   

            end
            
            
            trim = rectangle(:,center(1):end-(size(PSF,2)-center(1)));

         end
        function detector = CliptoDetector(order,wavelength,order_coeff,wave_coeff,cheby,p1,specorder)
        
         
            % convert each xypair into a set of 4 verticies for the pixel. Depends on PSF scaling etc..
            % from PSF convolution - know what one square "sample" is in mm. Does that still work after convolution?
            %% XY location pair for each point on wavelength vector
            xy_pair = zeros([size(order),2]);
            
            [~,beta] = max(order);
            center = mode(beta);
            
            % Convert lambda to X based on Zemax order fit
            % new x  =  interpolate (lam,x,lam')
            
%             xy_pair(center,:,1) = interp1( polyval(wave_coeff,((1:0.1:4096)-2048)/100)  , ((1:0.1:4096)-2048)/100   ,wavelength,'linear','extrap');
            
%             xy_pair(center,:,1) = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,wavelength,'linear','extrap');

            % Chebyshev 2D polynomial wavelength interpolation. 
             specorder = 40-specorder;
            if cheby == 1
                Inverse =1;
                nx = 7;
                nm = 7;
                ntotal = 39;
                npix = 4096;
                order0 = 113;
                finex=-25.48:0.01:25.48;
                fineorders = specorder*ones(size(finex)); % what spectral order are we at? 
                finechebs = Calculate_chebs(finex, fineorders+order0,Inverse,order0,ntotal,npix,nx,nm); 
                finewavelengths = JointFitFunc(p1,finechebs,fineorders+order0,nx,nm); 
                
                xy_pair(center,:,1) = interp1( finewavelengths,finex,wavelength,'linear','extrap');
                
            else
            % using 4th order, 1D polynomial needs interpolation. Inversion
            % is not possible. 
                 xy_pair(center,:,1) = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,wavelength,'linear','extrap');
            end
            
            % Inverting the 2nd order 1D polynomial 
%             a = wave_coeff(1);
%             b = wave_coeff(2);
%             c = wave_coeff(3);
% 
%             xy_pair(center,:,1) = (-b + sqrt(b^2-4*a*(c-wavelength)))/(2*a);
            
            
            % fprintf('Counts before detector: %.5f', sum(sum(order)))
            %% This makes the order flat
%              test = xy_pair(center,1,1).*ones(1,size(order,2));
%              xy_pair(center,:,2) = polyval(order_coeff,test);
            
            
           xy_pair(center,:,2) = polyval(order_coeff,xy_pair(center,:,1));
            
            % multiply by 100 to convert x,y mm to pixels then move 0,0 from center to bottom left
            xy_pair =xy_pair*100+4096/2;

            delta_w = diff(wavelength)/2;
            delta_w = [delta_w delta_w(end)];
                       
            edges_w = [wavelength(1)-delta_w(1) wavelength+delta_w];
            if cheby ==1
                edges_x = interp1( finewavelengths, finex ,edges_w,'linear','extrap');
            else
                edges_x = interp1( polyval(wave_coeff,-25.48:0.00001:25.48), (-25.48:0.00001:25.48) ,edges_w,'linear','extrap');
            end
%             edges_x = (-b + sqrt(b^2-4*a*(c-edges_w)))/(2*a);
            edges_x = edges_x*100+4096/2;
            
            
%             offset = xy_pair(center,2:size(xy_pair,2),1)-xy_pair(center,1:size(xy_pair,2)-1,1);            
%             offset = [offset offset(end)];
            
            yoffset = [center-1:-1:1 0 -(1:size(order,1)-center)];
            
            detector = zeros(4096+4,4096+4); % 4096x4096 + 2 extra for padding 
            
            %%% XY center pairs for entire order
            
            % xy_pair(vertical scalar index, horizontal vector of length, x and y)
            
            for jj = 1:size(xy_pair,2)
                for ii = [1:center-1 center+1:size(order,1)]
                    xy_pair(ii,jj,1) = xy_pair(center,jj,1);
%                     xy_pair(ii,jj,2) = xy_pair(center,jj,2)+ yoffset(ii)*offset(1);
                    xy_pair(ii,jj,2) = xy_pair(center,jj,2)+ yoffset(ii)*mean(diff(edges_x));

                end
            end
            
    
            for ii = [1:center-1 center center+1:size(order,1)]
                
                for jj =1:size(xy_pair,2)
%                     xlist{ii}{jj} = xy_pair(ii,jj,1) + [-1 1 1 -1 -1]*offset(jj)/2 ; %x starting at top left, clockwise
%                     ylist{ii}{jj} = xy_pair(ii,jj,2) + [1 1 -1 -1  1]*offset(1)/2; % y ***extra final point to draw box
                    xlist{ii}{jj} = [edges_x(jj) edges_x(jj+1) edges_x(jj+1) edges_x(jj) edges_x(jj)];
                    ylist{ii}{jj} = xy_pair(ii,jj,2) + [1 1 -1 -1  1]*mean(diff(edges_x))/2; % y ***extra final point to draw box
                    
                     % +1 to the x and y list to account for a 1 cell padd
                     % to remove edge clipping errors. 
                    [area,ind]=polyfillaa(2+xlist{ii}{jj},2+ylist{ii}{jj},4+4096, 4+4096);
                    ind = ind+1;
                    detector(ind) = detector(ind)+order(ii,jj).*area./sum(area);
                    clear ind area
                end
                
            end
            
            
            detector = detector(3:end-2,3:end-2); % remove padding
            detector = detector';
            
            
            
            % Plotting - be careful with this
%             figure(100);
%             hold on
%             
%             for jj = 1:200%size(order,2)
%             
%                 for ii = [1:center-1 center center+1:size(order,1)]
%             
%             
% %                     plot(xy_pair(:,jj,1),xy_pair(:,jj,2),'.k')
%             
%                     plot(xlist{ii}{jj},ylist{ii}{jj},'-','color',colors{rem(jj,length(colors)-1)+1})
%             
%             
%             
%                 end
%             end
            
        
        end
        function WriteFits(filename,image,headerinfo)
            
            import matlab.io.*
            DataType = 'single';
            % create fits file
            fptr = fits.createFile(filename);
            
            % Create Image extension
            fits.createImg(fptr,DataType,size(image));
            
            % Write image to file
            fits.writeImg(fptr,image);
            
            % Write additional keywords
            for ii = 1:size(headerinfo,1)
                fits.writeKey(fptr,headerinfo{ii,1},headerinfo{ii,2},headerinfo{ii,3});
                
            end    
            % close file
            fits.closeFile(fptr);
           
        end
        function [OutWavelength, OutSpectrum] = trimOrders (wavelength,spectrum,wave_coeff)
            
            ord = size(spectrum, 2); % what is the largert order (39)
            buffer_edge = 25; %set detector edge (buffered by 2mm)
            
            ind_lim1 = find(wavelength(:,ord) <= polyval(wave_coeff(ord,:,1),-buffer_edge),1,'last');
            ind_lim2 = find(wavelength(:,ord) >= polyval(wave_coeff(ord,:,1), buffer_edge),1,'first');
            max_length = (ind_lim2-ind_lim1); % find the length of the largest order with buffer

            for ii = 1:36

                ind1 = find(wavelength(:,ii) <= polyval(wave_coeff(ii,:,1), -buffer_edge),1,'last');
                ind2 = max_length + ind1; % ind2 is forced to match largest order size
%                 disp(ind2-ind1) %check all orders are the same size
                OutWavelength(:,ii) = wavelength(ind1:ind2,ii);
                OutSpectrum(:,ii) = spectrum(ind1:ind2,ii);
%                 disp(ii)
            end
            
        end
        function [y_grid] = resampleToGrid(x,y,x_grid)
            % this function cuts a new x vector and the max of the original and
            % extrapolates the short end of the new x vector to the short end of the
            % original x vector. The new y vector for use is generated by interpolating the
            % new_x,new_y onto the original x vector.
            
            extrap_scalar = min(y);
            ub = min(x(end,1),x_grid(end,1));
            ind = find(x<=ub,1,'last');
            x_safe = x(1:ind);
            y_safe = y(1:ind);
            
            y_grid = interp1(x_safe(:,1),y_safe(:,1),x_grid,'linear',extrap_scalar);
            
        end
    end
    
end

function [SR_interp]=SR_interpolation(StrehlSamples)

StrehlR = StrehlSamples(1,1);
StrehlI = StrehlSamples(2);
StrehlY = StrehlSamples(3);
StrehlJ = StrehlSamples(4);
StrehlH = StrehlSamples(5);
StrehlK = StrehlSamples(6);
StrehlL = StrehlSamples(7);
StrehlM = StrehlSamples(8);

% center of each band in wavelength space
% lets use mircrons
step = 1E-4;
Strehl = [StrehlR,StrehlI,StrehlY,StrehlJ,StrehlH,StrehlK,StrehlL,StrehlM];
Band_centers = [0.640,0.750,1.020,1.250,1.650,2.200,3.805,4.781];% from SOUL Excel document header information
xq = Band_centers(1,1):step:Band_centers(1,end);
SR = interp1(Band_centers, Strehl,xq);
SR_interp(:,1) = xq;
SR_interp(:,2) = SR;
end
function[x_cut,y_cut] = select_bandpass(x,y, bounds)
rng_u =max(bounds);
rng_l =min(bounds);
wavelength_limits = x<=rng_u & x>=rng_l; % Instrument band
y_cut = y(wavelength_limits);% cut the flux value at the band
x_cut = x(wavelength_limits);
area(x_cut,y_cut,'Facecolor','k','Linestyle','none')
line([x_cut(1,1) x_cut(1,1)], [0 y_cut(1,1)],'color','k')
line([x_cut(length(x_cut)) x_cut(length(x_cut))], [0 y_cut(length(y_cut))],'color','k')
alpha(0.15)
end
function [MagR,SRY,SRJ,SRR,SRI,SRH,SRK,SRL,SRM,fsHz,fluxphfrpix,bin,fluxphfr] = importStrehl...
  (filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [MAGR,SRR,SRI,SRY,SRJ,SRH,SRK,SRL,SRM,TOTALNM,BIN,MODAMPLD,FSHZ,FLUXPHFR,FLUXPHFRSA,FLUXPHFRPIX,MODES1,PHOTMODES,GAIN,TTGAIN]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [MAGR,SRR,SRI,SRY,SRJ,SRH,SRK,SRL,SRM,TOTALNM,BIN,MODAMPLD,FSHZ,FLUXPHFR,FLUXPHFRSA,FLUXPHFRPIX,MODES1,PHOTMODES,GAIN,TTGAIN]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [MagR,SRR,SRI,SRY,SRJ,SRH,SRK,SRL,SRM,Totalnm,bin,modamplD,fsHz,fluxphfr,fluxphfrsa,fluxphfrpix,modes1,Photmodes,gain,TTgain] = importfile('FLAO_Strehl.csv',2, 7);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/10/03 09:18:01

% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Allocate imported array to column variable names
MagR = dataArray{:, 1};
SRR = dataArray{:, 2};
SRI = dataArray{:, 3};
SRY = dataArray{:, 4};
SRJ = dataArray{:, 5};
SRH = dataArray{:, 6};
SRK = dataArray{:, 7};
SRL = dataArray{:, 8};
SRM = dataArray{:, 9};
Totalnm = dataArray{:, 10};
bin = dataArray{:, 11};
modamplD = dataArray{:, 12};
fsHz = dataArray{:, 13};
fluxphfr = dataArray{:, 14};
fluxphfrsa = dataArray{:, 15};
fluxphfrpix = dataArray{:, 16};
modes1 = dataArray{:, 17};
Photmodes = dataArray{:, 18};
gain = dataArray{:, 19};
TTgain = dataArray{:, 20};
end
function [tel] = read_orders(telnum,order_data)
% Read in raw Zemax spectral format data according to Macro written on 9/29/16
% Parse and Reshape data to produce upright 2D arrays corresponding to detector
% Add conversion in the future

%offsets measured by zooming into footprint diagram 12.2 11/1/17
%x offset = -10.578, y offset = -1.84

% order_data(:,2) = order_data(:,2)-10.578;
% order_data(:,3) = order_data(:,3)-1.84;  
% order_data(1:4:end,:)=[];

tel_w = reshape(order_data(telnum:3:end,1),[9,36])';
tel_y = reshape(-order_data(telnum:3:end,2),[9,36])';
tel_x = reshape(order_data(telnum:3:end,3),[9,36])';

% tel = zeros(39,9,3);

tel(:,:,1) = tel_x;
tel(:,:,2) = tel_y;
tel(:,:,3) = tel_w;

% Test figure
% figure(1)
% hold on
% plot(tel_x,tel_y,'ok')

end
function Sum=sumnd(Data)
%------------------------------------------------------------------------------
% sumnd function                                                     AstroStat
% Description: Return the global sum of a N-D matrix.
% Input  : - N-D matrix
% Output : - Value of global sum.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     March 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See also: meannd.m, minnd.m, maxnd.m, mediannd.m
% Example: Sum=sumnd(rand(4,4,4,4));
% Reliable: 1
%------------------------------------------------------------------------------

Sum = sum(Data(:));

%SizeData = size(Data);
%Vec = reshape(Data, [prod(SizeData), 1]);
%Sum = sum(Vec);
end 
function [areas,ret] = polyfillaa(px,py,sx,sy)



% Clip grid to the enclosing box

left = max(floor(min(px)),0);
right = min(floor(max(px)),sx-1);
bottom = max(floor(min(py)),0);
top = min(floor(max(py)),sy-1);
nx=right-left+1;
ny=top-bottom+1;


if nx < 1 || ny < 1
    
%     disp('clipping missed grid')
    ret = 0;
    areas = 0;
    return
end

%npix is the maximum possible number of clipped polys
npix=nx*ny;

if npix <= 0
%     disp('clipping missed grid')
    ret = 0;
    areas = 0;
    return
end

ind = 1;

areas = [];%zeros(1,npix);
ret = [];%zeros(1,npix);
% 
% figure(1)
% plot([px,px(1)],[py,py(1)],'-k')

for p = 1
    px1 = px;
    py1 = py;
    
    
    for j = bottom(p):top(p)
        for i =left(p):right(p)
            px_out = px1;
            py_out = py1;
%             polyclip(i,j,px_out,py_out)
            Polygon=[px_out;py_out]';
            unitPixel = [i,i+1,i+1,i;j,j,j+1,j+1]';
            
            clippedPolygon = sutherlandHodgman(Polygon,unitPixel);
%             figure(1)
%             hold on
                 
            
            if isempty(clippedPolygon) ==1
                continue 
                
            end
%             plot([clippedPolygon(:,1)',clippedPolygon(1,1)],[clippedPolygon(:,2)',clippedPolygon(1,2)],'-o')
           
            
            ret(ind)=i+j*sx;
            
            
            px_out = clippedPolygon(:,1);
            py_out = clippedPolygon(:,2);
            

                       
            areas(ind) = abs(sum(px_out.* circshift(py_out,-1)-py_out.*circshift(px_out,-1))./2);
            ind = ind+1;
        end
    end
end


zero_ind = find(areas==0);
areas(zero_ind)=[];
ret(zero_ind)=[];

end
function clippedPolygon = sutherlandHodgman(subjectPolygon,clipPolygon)
%The inputs are a table of x-y pairs for the verticies of the subject
%polygon and boundary polygon. (x values in column 1 and y values in column
%2) The output is a table of x-y pairs for the clipped version of the 
%subject polygon. 
%% Helper Functions
 
    %computerIntersection() assumes the two lines intersect
    function intersection = computeIntersection(line1,line2)
 
        %this is an implementation of
        %http://en.wikipedia.org/wiki/Line-line_intersection
 
        intersection = zeros(1,2);
 
        detL1 = det(line1);
        detL2 = det(line2);
 
        detL1x = det([line1(:,1),[1;1]]);
        detL1y = det([line1(:,2),[1;1]]);
 
        detL2x = det([line2(:,1),[1;1]]);
        detL2y = det([line2(:,2),[1;1]]);
 
        denominator = det([detL1x detL1y;detL2x detL2y]);
 
        intersection(1) = det([detL1 detL1x;detL2 detL2x]) / denominator;
        intersection(2) = det([detL1 detL1y;detL2 detL2y]) / denominator;
 
    end %computeIntersection
 
    %inside() assumes the boundary is oriented counter-clockwise
    function in = inside(point,boundary)
 
        pointPositionVector = [diff([point;boundary(1,:)]) 0];
        boundaryVector = [diff(boundary) 0];
        crossVector = cross(pointPositionVector,boundaryVector);
 
        if ( crossVector(3) <= 0 )
            in = true;
        else
            in = false;
        end
 
    end %inside
 
%% Sutherland-Hodgman Algorithm
 
    clippedPolygon = subjectPolygon;
    numVerticies = size(clipPolygon,1);
    clipVertexPrevious = clipPolygon(end,:);
 
    for clipVertex = (1:numVerticies)
 
        clipBoundary = [clipPolygon(clipVertex,:) ; clipVertexPrevious];
 
        inputList = clippedPolygon;
 
        clippedPolygon = [];
        if ~isempty(inputList)
            previousVertex = inputList(end,:);
        end
 
        for subjectVertex = (1:size(inputList,1))
 
            if ( inside(inputList(subjectVertex,:),clipBoundary) )
 
                if( not(inside(previousVertex,clipBoundary)) )  
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);
                end
 
                clippedPolygon(end+1,1:2) = inputList(subjectVertex,:);
 
            elseif( inside(previousVertex,clipBoundary) )
                    subjectLineSegment = [previousVertex;inputList(subjectVertex,:)];
                    clippedPolygon(end+1,1:2) = computeIntersection(clipBoundary,subjectLineSegment);                            
            end
 
            previousVertex = inputList(subjectVertex,:);
            clipVertexPrevious = clipPolygon(clipVertex,:);
 
        end %for subject verticies                
    end %for boundary verticies
end %sutherlandHodgman
function F=circ_gauss(X,Y,Sigma,center)
%--------------------------------------------------------------------------
% circ_gauss function                                                General
% Description: Calculate 2D circular Gaussian in a 2-D grid.
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D Gaussian.
%          - same as the x-ccordinates, but for the y-axis.
%          - Sigma of the Gaussian or [SigmaX, SigmaY] in case sigma
%            is different for each axis.
%            By default SigmaY=SigmaX.
%            If empty matrix use default.Si
%          - Center of the Gaussian [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - Maximum radius of Gaussian behond to set it to zero.
%            Default is Inf.
%            MaxRad is measured from the center of the kernel and not
%            the center of the Gaussian.


% Example: 
%          F=circ_gauss(MatX,MatY,[1],[0 0]);
%          surface(F);
%--------------------------------------------------------------------------



SigmaX = Sigma(1);
SigmaY = Sigma(2);

X0 = center(1);
Y0 = center(2);

F = 1./(2.*pi.*SigmaX.*SigmaY) .*exp(-1./(2.).* ((X-X0).^2./SigmaX.^2 +(Y-Y0).^2./SigmaY.^2));

F = F./sum(sum(F));

% set elements outside MaxRad to zero:
% if (~isinf(cutoff)),
%    MatR = sqrt(X.^2 + Y.^2);
%    I = find(MatR>cutoff);
%    F(I) = 0;
% end
% 
% if (isnan(Norm)),
%    % do not normalize
% else
%    F = Norm.*F./sumnd(F);
% end
end