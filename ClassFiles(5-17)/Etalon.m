classdef Etalon < Spectra
    properties
        FSR
        Finesse
        Length
        R
        refpeak
        FSR_measured
        FSR_estimated
        RV
    end
    
    methods
        function [obj] = Etalon(EtalonType,FSR,scale,TestRV)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
            else
                % need to put conditions for user specified inputs here
            end
            
            obj = LoadEtalon(obj,EtalonType,FSR,scale); %Andrew added this 5/23/17
            obj.RV = TestRV;
            obj.Wavelength = obj.DopplerShift;
        end
        
        
        function [obj] = LoadFP(obj)
            global calibfile
            load(calibfile); %M Foleys arbitrary lorentzians
            obj.Wavelength = (Lorentzian(1,:)*1e10)'; % Angstroms 
            obj.Counts = (Lorentzian(2,:)*0.5e7/1800)'; % BS scale factor
        end       
        function [obj] = LoadEtalon(obj,type,FSR,scale)
            R1 = 0.93; % Mirror 1 reflectivity
            R2 = 0.93; % Mirror 2 reflectivity
%             l = 14.3e-3;  % Cavity length (m)
%             FSR = 5; %GHz
            % gauss params.
            sep = 1; %nm separation
            sigma = 0.02*1.5/2.35; % sigma nm;
            
            if strcmp(type,'gauss')               
                
            [wavelength,counts] = GaussianSpec(sep,sigma);
            obj.Wavelength(:,1) = wavelength*10; % Angstroms
            obj.Counts(:,1) = counts; % BS scale factor            
            
            elseif strcmp(type,'delta')               
                
            [wavelength,counts] = DeltaSpec(sep);
            obj.Wavelength(:,1) = wavelength*10; % Angstroms
            obj.Counts(:,1) = counts; % BS scale factor  
                
            
            
            else 
            
            [wavelength,Intesnity,p,cent] = GenerateEtalon (R1,R2,FSR,type,scale);
            
            %-----Assign values to Etalon object-----%
            obj.Wavelength(:,1) = wavelength*1e10; % Angstroms
            obj.Counts(:,1) = Intesnity; % BS scale factor
            obj.Length = p(1,1); %cavity length
            obj.R = p(1,2); % R1 reflectivity
            obj.Finesse = p(1,3); % Finesse 
            obj.FSR = p(1,4); %FSR in Ghz
%             obj.dlambda = []; %centers in wavelength space
            obj.refpeak = cent;
            end
        end    
        function [] = PlotEtalon (obj)
            close all
%             figure()
%             plot(obj.Wavelength/10,obj.dlambda)
                        
            figure()
            plot(obj.Wavelength/10, obj.Counts , 'linewidth' , 0.5 , 'color' , [0.8 , 0.1 ,0.8]);
            ylim([0 1.2])
            
% %             window = 5*0.1*1E9;
%             center = mean(obj.Wavelength)/10;
%             lb = center-window;
%             ub = center+window;
% %             xlim([lb ub])
            
            disp(['FSR = ',num2str(obj.FSR),'GhZ'])
            disp(['F = ',num2str(obj.Finesse)])
        end             
        function [obj] = FSRCalc (obj)
            
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
       function [Counts3] = Scale2Star (Counts1, Counts2)
       
           A = max(max(Counts1)); %max of the Etalon spectrum;
           B = max(max(Counts2)); %max of the Stellar spectrum
           ScaleFactor = B/A; %What is the ratio of the stat to cal
           Counts3 = Counts1*ScaleFactor*10; %Scale the calibration source to the star
            
       end
        
    end
end


function [lambda,I,p,cent] = GenerateEtalon (R1,R2,FSR,type,scale)
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
l = 1/(FSR*1E9*(2*n/c));
theta = 0; %normal incidence

%-----Calculations-----%
delta = (2*pi./lambda)*2*n*l*cos(theta)+2*phi; % ck/FSR
f_coeff = 4*sqrt(R1*R2)/((1-R1)*(1-R2));
Finesse = pi*(R1*R2)^(1/4)/(1-sqrt(R1*R2));
% f_coeff_check = (2*finesse/pi).^2;
I = 1 ./ (1 + f_coeff * sin(delta/2).^2);

% dlambda = (lambda.^2)./(2*l*n+lambda);
% lambda_spacing = 1E9*mean(dlambda);


if strcmp(type, 'science')==1
    cent = 0;
elseif strcmp(type, 'ref')==1
    %Reference etalon
    lb = 975; 
    ub = lb+0.02;
    lambda_ref = (lb:0.000001:ub)*1e-9;
    delta_ref = (2*pi*2*n*l./(lambda_ref)+2*phi);
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
        delta_ref = (2*pi*2*n*l./(lambda_ref)+2*phi);
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

%-----Outputs-----%
p = [l,R1,Finesse,FSR]; %parameters
end

function [wave, spec] = GaussianSpec(sep,width)


wave = 920:0.0001:1400; 
spec =zeros(size(wave));

N = (wave(end)-wave(1))/sep;



for ii = 1:N
   
    center(ii) = wave(1)+(ii)*sep;
       
    
    spec = spec + normpdf(wave,center(ii),width);
    
    
end
end
function [wave, spec] = DeltaSpec(sep)

wave = 920:5e-4:1400; 
spec =zeros(size(wave));

pwave = 920:sep:1400; 

for ii = 1:length(pwave)
    ind = find(wave==pwave(ii));
    spec(ind)=1;

end

end

