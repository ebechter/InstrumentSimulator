classdef Imager < Instrument
    properties
        psf %fwhm
        pixelPitch %pixel size
        detectorDimensions % format in pixels
        plateScale %mas/micron or as/mm
        apDiameter
        blockFrac
        frame % sample frame
        q % FWHM sampling
    end
    
    methods
        function [obj] = Imager(type)
            %% Pre Initialization %%
            %% NIR II = front and back AR coating.
            %% NIR  = single surface NIR coating.
            blockFrac = [];
            apDiameter = [];
            name = [];
            plateScale = [];
            pixelPitch = [];
            detectorDimensions = [];
            bandPass = [];
            psf= [];
            % Any code not using output argument (obj)
            if nargin == 0 || strcmp(type,'ANDOR') == 1
                %Need to define an optical model structure and rules.
                
                
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{4} = struct('name','L1_1','type','triplet lens','coatingName','RMI_BASF51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257,'internalTrans',0.998*0.996*0.999);
                opticalModel{5} = struct('name','L1_2','type','triplet lens','coatingName','RMI_FPL51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257);
                opticalModel{6} = struct('name','ADC1_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{7} = struct('name','ADC1_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{8} = struct('name','ADC2_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{9} = struct('name','ADC2_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{10} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{11} = struct('name','Shortpass','type','Dichroic','coatingName','AlluxaSP_T','number',1,'angle','45','efficiency',[],'surfaceQuality',1/10,'focalLength',[]);
                opticalModel{12} = struct('name','L4_2','type','triplet lens','coatingName','RMI_FPL51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257);
                opticalModel{13} = struct('name','L4_1','type','triplet lens','coatingName','RMI_BASF51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257,'internalTrans',0.998*0.996*0.999);
                opticalModel{14}= struct('name','ANDOR','type','detector','coatingName','AndorZyla42','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                
                name = 'Andor Channel';
                plateScale = [];
                pixelPitch = 6.5e-6;
                detectorDimensions = [2048,2048];
                bandPass = [800,1000];
                psf= 47e-6;
                
                
            elseif strcmp(type,'WFC') == 1
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{4} = struct('name','WFCL1','type','doublet lens','coatingName','NIRII','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/4,'focalLength',250);
                %opticalModel{5} = struct('name','Basler','type','detector','coatingName','Basler','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                
                name = 'Wide Field Channel';
                bandPass = [800,1100];
                pixelPitch = 5.5e-6;
                detectorDimensions = [2048,2048];
                plateScale =[];
                psf = 12.5e-6; %950nm
                
                
            elseif strcmp(type,'FIBER') == 1

                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{4} = struct('name','L1_1','type','triplet lens','coatingName','RMI_BASF51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257,'internalTrans',0.998*0.996*0.999);
                opticalModel{5} = struct('name','L1_2','type','triplet lens','coatingName','RMI_FPL51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257);
                opticalModel{6} = struct('name','ADC1_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{7} = struct('name','ADC1_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{8} = struct('name','ADC2_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{9} = struct('name','ADC2_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{10} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{11} = struct('name','Shortpass','type','Dichroic','coatingName','AlluxaSP_R','number',1,'angle','45','efficiency',[],'surfaceQuality',1/10,'focalLength',[]);
                opticalModel{12} = struct('name','Longpass','type','Dichroic','coatingName','AlluxaLP_R','number',1,'angle','45','efficiency',[],'surfaceQuality',1/10,'focalLength',[]);
                opticalModel{13} = struct('name','L2_1','type','triplet lens','coatingName','RMI_BASF51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257,'internalTrans',0.998*0.996*0.999);
                opticalModel{14} = struct('name','L2_2','type','triplet lens','coatingName','RMI_FPL51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257);
                opticalModel{15} = struct('name','L3_1','type','triplet lens','coatingName','ICOS_LAL14','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',20,'internalTrans',0.998*0.97*0.996);
                opticalModel{16} = struct('name','L3_2','type','triplet lens','coatingName','ICOS_BSL7','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',20);
                

                name = 'Fiber Channel';
                bandPass = [965,1300];
% %                 bandPass = [600,1500];
                pixelPitch = 0.5e-6;
                detectorDimensions = [100,98];
                plateScale = [];
                psf = 3.6e-6;% 980nm;
                
            elseif strcmp(type,'QUADCELL') == 1
                
               opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{4} = struct('name','L1_1','type','triplet lens','coatingName','RMI_BASF51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257,'internalTrans',0.998*0.996*0.999);
                opticalModel{5} = struct('name','L1_2','type','triplet lens','coatingName','RMI_FPL51','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/10,'focalLength',257);
                opticalModel{6} = struct('name','ADC1_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{7} = struct('name','ADC1_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{8} = struct('name','ADC2_1','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[],'internalTrans',0.998*0.9965*0.997);
                opticalModel{9} = struct('name','ADC2_2','type','triplet prism','coatingName','ICOS_ADC','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{10} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[],'surfaceQuality',1/20,'focalLength',[]);
                opticalModel{11} = struct('name','Shortpass','type','Dichroic','coatingName','AlluxaSP_R','number',1,'angle','45','efficiency',[],'surfaceQuality',1/10,'focalLength',[]);
                opticalModel{12} = struct('name','Longpass','type','Dichroic','coatingName','AlluxaLP_T','number',1,'angle','45','efficiency',[],'surfaceQuality',1/10,'focalLength',[]);
                opticalModel{13}= struct('name','L5','type','doublet lens','coatingName','NIRII','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/4,'focalLength',[]);
                %opticalModel{11} = struct('name','QuadFilter','type','bandpass filter','coatingName','Imagefilter','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                opticalModel{14}= struct('name','QuadCell','type','detector','coatingName','NIRQuadCell','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                
                name = 'Quad Cell Channel';
                bandPass = [1300,1600];
                pixelPitch = [0.25e-3];
                detectorDimensions = [2,2];
                plateScale = [];
                psf = 51e-6; %FWHM at 1550nm
                
            elseif strcmp(type,'LBT') == 1
                opticalModel{1} = struct('name','Primary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
                opticalModel{2} = struct('name','Secondary','type','mirror','coatingName','LBT_T_Al','number',1,'angle',[],'efficiency',[]);
                opticalModel{3} = struct('name','Tertiary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [400,2000];
                plateScale = [];
                blockFrac = 0.108; % value listed as 0.889
                apDiameter = 8.22; % meters (IR effective collecting area, supposedly)
                name = 'LBT';
                pixelPitch = 6.5e-6;
                detectorDimensions = [512,512];
                psf = 1.7112e-05; %scaled using ANDOR camera
                
                
            elseif strcmp(type,'LBTI') == 1
                opticalModel{1} = struct('name','Entrance Window','type','dichroic','coatingName','AlluxaEntWT','number',1,'angle','15','efficiency',[]);
%                 opticalModel{1} = struct('name','Entrance Window','type','dichroic','coatingName','LBTI_ZnSeT','number',1,'angle',[],'efficiency',[]);
                opticalModel{2} = struct('name','Ellipse Mirror','type','elliptical mirror','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','Pupil mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','Roof mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{5} = struct('name','Switch Box mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{6} = struct('name','Exit window','type','mirror','coatingName','RMI_exit','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [400,2000];
                pixelPitch = 6.5e-6;
                detectorDimensions = [512,512];
                plateScale = [];
                name = 'LBTI';
                psf = 1.7112e-05; %scaled using ANDOR camera
                
            elseif strcmp(type,'Filter') == 1
                opticalModel{1} = struct('name','Filter','type','bandpass','coatingName','FB1330_12','number',1,'angle',[],'efficiency',[]);
                bandPass = [600,1500];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                name = 'Filter';
                
                
            elseif strcmp(type,'flat') == 1
                opticalModel{1} = struct('name','ones','type','none','coatingName','ones','number',1,'angle',[],'efficiency',[]);
                bandPass = [970,1270];
                plateScale = [];
                name = 'flat';
                pixelPitch = 6.5e-6;
                detectorDimensions = [512,512];
                psf = 1.7112e-05; %scaled using ANDOR camera
            
            elseif strcmp(type,'demonstrator') == 1
                opticalModel{1} = struct('name','demo','type','none','coatingName','demo','number',1,'angle',[],'efficiency',[]);
                opticalModel{2} = struct('name','Input power','type','thor lens','coatingName','NIRII','number',1,'angle',[],'efficiency',[],'surfaceQuality',1/4);
                opticalModel{3}= struct('name','QuadCell','type','detector','coatingName','NIRQuadCell','number',1,'angle',[],'efficiency',[],'surfaceQuality',[],'focalLength',[]);
                
                bandPass = [970,1065];
                plateScale = [];
                name = 'demonstrator';
                pixelPitch = 6.5e-6;
                detectorDimensions = [512,512];
                psf = 1.7112e-05; %scaled using ANDOR camera
                
            end  
            
            
            current_path = pwd;
            if strcmp(current_path(2:8),'Volumes')==1
               curveDirectory = '/Volumes/Software/Simulator/RefFiles/Curves/Instrument/';

            elseif strcmp(current_path(2:4),'afs')==1
                
                curveDirectory = '/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/RefFiles/Curves/Instrument/';
            else
                
                curveDirectory = [current_path(1:2) '\Simulator\RefFiles\Curves\Instrument\'];
            end
            
            obj.bandPass = bandPass;
            obj.pixelPitch = pixelPitch;
            obj.detectorDimensions = detectorDimensions ;
            obj.opticalModel = opticalModel;
            obj.plateScale = plateScale;
            obj.psf = psf;
            obj.blockFrac = blockFrac;
            obj.apDiameter = apDiameter;
            obj.name = name;
            obj.q = psf./pixelPitch;
            obj.frame = obj.genPSF(bandPass.*1e-9,obj.q,detectorDimensions); % meters, q, pixels
            
            [obj] = loadOpticalModelCurves(obj,curveDirectory);
            [obj] = trimThroughput(obj);
            [obj] = intThroughput(obj);
        end
    end
    
    methods(Static)
        
        function [airy] = genPSF(bandPass, q, npup)
        
            lambda = mean(bandPass); % in meters
            
            alpha = 0.11;% blocking parameters
            
            npup = min(npup);
            
            [airy] = diffraction_pattern(q,lambda,npup,alpha); % compute diffraction pattern
            
        end
        
    end
end


function C = circ(N,R)
C = zeros(N);

[X,Y] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);

X = X+0.5; Y=Y+0.5;

Z = sqrt(X.^2 + Y.^2);

C(Z<R) = 1;

end
function [airy] = diffraction_pattern(q,lambda,npup,alpha)
            
            %%-------------------
            % upscale computation
            %%-------------------
            
            q = 2*q; % gives higher sampling
            
            npup = 2*npup; % keeps grid size same after binning by factor of 2 (factor of 2 is used to generate upscale simluated nyquist PSFs by minimum amount) 
            
            %%-----------------------
            % calculate airy pattern
            %%-----------------------
            
            pup_out = circ(npup,npup/q/2);
            
            pup_in = circ(npup,alpha*npup/q/2);
            
            pupil = pup_out-pup_in;
            
            opd = ones(size(pupil));
            
            P = pupil .* exp(1i * 2*pi/lambda * opd);
            
            psf = fftshift(fft2(P));
            
            psf = psf .* conj(psf);
            
            airy = psf ./ sum(psf(:));
            
            
            %%----------
            % Bin down
            %%----------
            
            p = 2; % default by 2
            
            q = 2; % default by 2
            
            [m,n]=size(airy); %M is the original matrix

            airy=sum(reshape(airy,p,[]) ,1 );
            
            airy=reshape(airy,m/p,[]).'; %Note transpose

            airy=sum(reshape(airy,q,[]) ,1);
            
            airy=reshape(airy,n/q,[]).'; %Note transpose

        end