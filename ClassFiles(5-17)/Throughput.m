classdef Throughput
    properties
        LBT
        EntWindow
        LBTI
        WFC
        CommonCh
        FiberCh
        ImageCh
        QuadCh
        Spec
        Grating
        FiberLink
        FullBand
        LBTIAlignmentCh
        CombinedFPCh
    end
    methods
        
        function[obj] = LBTMirrors(obj)
            global curve_dir
            %% ----------Load Curves----------%
            filename = 'mirr_PS_LBT.mat';
            file = strcat(curve_dir,filename);
            load(file)
            filename = 'mirr_T_LBT.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Central Obstruction losses----------%
            LBT_M1_eff = mirr_PS_LBT(:,2); %coating reflectivity
            LBT_M2_eff = mirr_PS_LBT(:,2);
            LBT_M3_eff = mirr_T_LBT(:,2);
            M1_radius = 8.22/2;
            M2_radius = 0.9124/2;
            A2=pi*(M2_radius)^2; %area of M2
            A1=pi*(M1_radius)^2; %area of M1
            Cen_O = 1-A2/A1; % Central obstruction ~ 11%
            %% ----------Multiply Curves----------%
            LBT_eff = LBT_M1_eff.*LBT_M2_eff.*LBT_M3_eff*Cen_O;
            %% ----------Assign Object properties----------%
            obj.LBT(:,1) = mirr_PS_LBT(:,1);
            obj.LBT(:,2) = LBT_eff;
        end
        function[obj] = LBTIRetro(obj)
            global curve_dir
            %% ----------Load Curves----------%
            filename = 'mirr_PS_LBT.mat';
            file = strcat(curve_dir,filename);
            load(file)
            filename = 'mirr_T_LBT.mat';
            file = strcat(curve_dir,filename);
            load(file)
            Retro = 0.92;
            %% ----------Central Obstruction losses----------%
            LBT_M1_eff = mirr_PS_LBT(:,2); %coating reflectivity
            LBT_M2_eff = mirr_PS_LBT(:,2);
            LBT_M3_eff = mirr_T_LBT(:,2);
            M1_radius = 8.22/2;
            M2_radius = 0.9124/2;
            A2=pi*(M2_radius)^2; %area of M2
            A1=pi*(M1_radius)^2; %area of M1
            Cen_O = 1-A2/A1; % Central obstruction ~ 11%
            %% ----------Multiply Curves----------%
            LBTI_eff = LBT_M3_eff.*LBT_M2_eff.*Retro.*LBT_M2_eff.*LBT_M3_eff*0.5*0.95;
            %% ----------Assign Object properties----------%
            obj.LBTIAlignmentCh(:,1) = mirr_PS_LBT(:,1);
            obj.LBTIAlignmentCh(:,2) = LBTI_eff;
        end
        function[obj] = EntranceWindow(obj,entrance_dichroic)
            
            global curve_dir
            
            if strcmp(entrance_dichroic,'ECI') ==1 
            %% ----------Load Curves----------%
            filename = 'EntWTrans';
            file = strcat(curve_dir,filename);
            load(file)
            LBTI_bs_wave = EntWTrans(:,1);
            LBTI_bs_trans = EntWTrans(:,2);
            
            elseif strcmp(entrance_dichroic,'sim') ==1 
            %% ----------Generate Curves----------%
            type = 'LBTI';
            [LBTI_bs_wave,LBTI_bs_trans]=Science_splitter(type);
            
            elseif strcmp(entrance_dichroic,'existing') ==1 
            %% ----------Load Curves----------%
            filename = 'LBTI_existing.mat';
            file = strcat(curve_dir,filename);
            load(file)
            xq = 0.5:0.001:2;
            vq=interp1(LBTI_existing(:,1),LBTI_existing(:,3),xq);
            LBTI_bs_wave = 1000*xq;
            LBTI_bs_trans = vq/100;
            end
            
            %% ----------Assign Object properties----------%
            obj.EntWindow(:,1)=LBTI_bs_wave;
            obj.EntWindow(:,2)=LBTI_bs_trans;
        end
        function[obj] = EntranceWindowWFS(obj,entrance_dichroic)
            global curve_dir
            
            if strcmp(entrance_dichroic,'ECI') ==1 
            %% ----------Load Curves----------%
            filename = 'EntWRefl.mat';
            file = strcat(curve_dir,filename);
            load(file)
            LBTI_bs_wave = EntWRefl(:,1);
            LBTI_bs_trans = EntWRefl(:,2);
            
            elseif strcmp(entrance_dichroic,'sim') ==1 
            %% ----------Generate Curves----------%
            type = 'AO';
            [LBTI_bs_wave,LBTI_bs_trans]=Science_splitter(type);
            
            elseif strcmp(entrance_dichroic,'existing') ==1
            %% ----------Load Curves----------%
            filename = 'LBTI_existing.mat';
            file = strcat(curve_dir,filename);
            load(file)
            xq = 0.5:0.001:2;
            vq=interp1(LBTI_existing(:,1),LBTI_existing(:,2),xq);
            LBTI_bs_wave = 1000*xq;
            LBTI_bs_trans = vq/100;
            end
            
            %% ----------Assign Object properties----------%
            obj.EntWindowWFS(:,1)=LBTI_bs_wave;
            obj.EntWindowWFS(:,2)=LBTI_bs_trans;
        end
        function[obj] = WFSOptics(obj,AO)
            global curve_dir
            Optics = 0.9635^11;
            
            %% ----------Load Curves----------%
            filename = 'detOCAM.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'detCCD39.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Multiply Curves----------%
            FLAO = Optics.* det_CCD39(:,2);
            SOUL = Optics.* det_OCAM(:,2);
            %% ----------Assign Object properties----------%
            if nargin > 1
                if AO == 1
                    obj.WFS(:,1)= det_CCD39(:,1);
                    obj.WFS(:,2)= FLAO;
                elseif AO ==2 || AO ==3 || AO ==4
                    obj.WFS(:,1)= det_OCAM(:,1);
                    obj.WFS(:,2)= SOUL;
                end
            else
                obj.WFS(:,1)= det_OCAM(:,1);
                obj.WFS(:,2)= SOUL;
                disp('SOUL WFS QE is used by default')
            end
            
            
            
        end
        function[obj] = LBTIOptics(obj)
            global curve_dir
            %% ----------Load Curves----------%
            
            filename = 'mirrGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'windSil.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'windSil.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            %% ----------Multiply Curves----------%
            [wav,Tput]=multiply_curves(mirr_Gold(:,1),mirr_Gold(:,2).^4,wind_Sil(:,1),wind_Sil(:,2));
            %% ----------Assign Object properties----------%
            obj.LBTI(:,1) = wav;
            obj.LBTI(:,2) = Tput;
            
        end
        function[obj] = CommonChOptics(obj)
            global curve_dir
            %% ----------Load Curves----------%
            filename = 'lensNIR.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'mirrGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(mirr_Gold(:,1),mirr_Gold(:,2).^4,lens_NIR(:,1),lens_NIR(:,2).^3);
            %% ----------Assign Object properties----------%
            obj.CommonCh(:,1) = wav;
            obj.CommonCh(:,2) = Tput;
        end
        function[obj] = FiberChOptics(obj)
            global curve_dir
            %% ----------Load lens Curves----------%           
            filename = 'lensNIR.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            %% ----------Load dichroic Curves----------%
            filename = 'FiberDRefl';
            file = strcat(curve_dir,filename);
            load(file)
            fiber_bs_wave = FiberDRefl(:,1);
            fiber_bs_trans = FiberDRefl(:,2);
            
            filename = 'QuadDRefl';
            file = strcat(curve_dir,filename);
            load(file)
            quad_bs_wave = QuadDRefl(:,1);
            quad_bs_trans = QuadDRefl(:,2);

            %% ----------Generate Curves----------%
            %This is no longer used now that we have quotes
            %type = 'fiber';
            %[sci_bs_wave,sci_bs_trans]=Science_splitter(type);
            
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(fiber_bs_wave,fiber_bs_trans,quad_bs_wave,quad_bs_trans);
            [wav,Tput]= multiply_curves(wav,Tput,lens_NIR(:,1),lens_NIR(:,2).^2);
            
            %% ----------Assign Object properties----------%
            obj.FiberCh(:,1) = wav;
            obj.FiberCh(:,2) = Tput;
        end
        function[obj] = ImageChOptics(obj)
            global curve_dir
            %% ----------Load Curves----------%
            filename = 'lensNIR.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'mirrGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'detZyla42.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Generate Curves----------%
            type = 'guide';
            [sci_bs_wave,sci_bs_trans]=Science_splitter(type);
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(sci_bs_wave,sci_bs_trans,lens_NIR(:,1),lens_NIR(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,mirr_Gold(:,1),mirr_Gold(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,det_Zyla42(:,1),det_Zyla42(:,2));
            %% ----------Assign Object properties----------%
            obj.ImageCh(:,1) = wav;
            obj.ImageCh(:,2) = Tput;
        end
        function[obj] = WFCOptics(obj)
            global curve_dir
            %% ----------Load Curves----------%
            filename = 'lensNIR.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'mirrGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'detZyla55.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(mirr_Gold(:,1),mirr_Gold(:,2),lens_NIR(:,1),lens_NIR(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,det_Zyla42(:,1),det_Zyla42(:,2));
            %% ----------Assign Object properties----------%
            obj.WFC(:,1) = wav;
            obj.WFC(:,2) = Tput;
        end
        function[obj] = QuadOptics(obj)
            global curve_dir
            %% ----------Generate Curves----------%
            type = 'guide';
            [sci_bs_wave,sci_bs_trans]=Science_splitter(type);
            %% ----------Load Curves----------%
            filename = 'lensNIR.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'mirrGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'detQuad.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'bsSpass.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(sci_bs_wave,sci_bs_trans,bs_Spass(:,1),bs_Spass(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,mirr_Gold(:,1),mirr_Gold(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,lens_NIR(:,1),lens_NIR(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,det_Quad(:,1),det_Quad(:,2));
            %% ----------Assign Object properties----------%
            obj.QuadCh(:,1) = wav;
            obj.QuadCh(:,2) = Tput;
        end
        function[obj] = SpecOptics(obj)
            global curve_dir
            %% ----------Load Curves----------%
            %fathom gold file format: 
            % WL(nm),10deg r-pol,10deg s-pol,10deg p-pol,25deg r-pol,25deg s-pol,25deg p-pol
            filename = 'FathomGold.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            Spol = (FathomGold(:,3)/100).^6.*(FathomGold(:,6)/100).^5;
            Ppol = (FathomGold(:,4)/100).^6.*(FathomGold(:,7)/100).^5;
            Rpol = (Spol+Ppol)/2;
            
            filename = 'detH4RG.mat';
            file = strcat(curve_dir,filename);
            load(file)
            p1 = polyfit(det_H4RG(:,1),det_H4RG(:,2),10);
            H4RG = polyval(p1,det_H4RG(:,1)); % smoothed out crappy H4RG QE data
            
            filename = 'grt_xdisp.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            %spectral filter
            specfilter = 0.95;
            
            %% ----------Multiply Curves----------%
            [wav,Tput]= multiply_curves(FathomGold(:,1),Rpol,grt_xdisp(:,1),grt_xdisp(:,2));
            [wav,Tput]= multiply_curves(wav,Tput,det_H4RG(:,1),H4RG);
            %% ----------Assign Object Properties----------%
            obj.Spec(:,1) = wav;
            obj.Spec(:,2) = Tput.*specfilter;
%             %% ----------Old Spectrograph Load Curves----------%
%             filename = 'mirrGold.mat';
%             file = strcat(curve_dir,filename);
%             load(file)
%             
%             filename = 'detH4RG.mat';
%             file = strcat(curve_dir,filename);
%             load(file)
%             
%             filename = 'grt_xdisp.mat';
%             file = strcat(curve_dir,filename);
%             load(file)
%             %% ----------Multiply Curves----------%
%             [wav,Tput]= multiply_curves(mirr_Gold(:,1),mirr_Gold(:,2).^11,grt_xdisp(:,1),grt_xdisp(:,2));
%             [wav,Tput]= multiply_curves(wav,Tput,det_H4RG(:,1),det_H4RG(:,2));
%             %% ----------Assign Object Properties----------%
%             obj.Spec(:,1) = wav;
%             obj.Spec(:,2) = Tput;
            
        end
        function[obj] = R6Grating(obj)
            %% ----------Load Curves----------%
            global gratingfile
            load(gratingfile)
            %% ----------Assign Object Properties----------%
            if isempty (obj.Polarization  ==1) || obj.Polarization(1,3) == 0
            obj.Grating=GratingEff_new;
            else
               pfrac = obj.Polarization(1,2);
               sfrac = 1-pfrac;
               dop =  obj.Polarization(1,1);
               
               %offset s and p in wavelength space
               offset = 0.5; %1/4 to 1/2 nm looks about right from measured data 
               peff = GratingEff_new(:,2:40);
               wave = GratingEff_new(:,1);
               wave_s = wave + offset;
               seff = interp1(wave_s,peff,wave);
               
               %rescale amplitudes
               scale = 1.135;
%              scale = 1.5;
               peff = scale*peff;
               seff = (2-scale)*seff;              
               
               unpolarized = 0.5*(seff+peff);
               polarized = dop*(seff*sfrac+peff*pfrac)+unpolarized*(1-dop);
               obj.Grating(:,1)=GratingEff_new(:,1);
               obj.Grating(:,2:40)=polarized;
               
            end
            
        end
        function[obj] = iLocFiberLink(obj,coupling,lambda)
            %----------Light Losses realted to Fiber injection parameters----------%
            Fiber_Fresnel = 0.92; %fresnel refelction at an uncoated fiber tip amd exit
            
            OverlapIntegral = interp1(1000*lambda,coupling,obj.FullBand','linear','extrap'); %mode mismatch overlap integral
            
            %----------Strehl Ratio Terms----------%
            %Dynamic Strehl ratio is the AO delivered strehl ratio
%             SR_AO = interp1(1000*SR_interp(:,1),SR_interp(:,2)/100,wav_f);%Strehl ratio resampled onto the fiber wavelength vector
            
            %Static Strehl ratio delivered by the optics
            [SR_Static]=Wagner();
            
            %----------Light Losses realted to Fiber propogation----------%
            Avim_connector = 0.98;
            number = 3;
            Avim = Avim_connector^number;
            Spec_Fiberlink = Avim*0.99;
            FWR_Fiberlink = Avim*0.01;
            dBLoss = 1-(10^(1.2/10)*45/1000); % fibercore spec sheet
            OverlapMod = (SR_Static.*OverlapIntegral)';
            
            %----------Assign Object values----------%
            obj.FiberLink(:,1) = obj.FullBand';
            obj.FiberLink(:,2) = Fiber_Fresnel*Spec_Fiberlink.*OverlapMod.*dBLoss;
            obj.FiberLink(:,3) = Fiber_Fresnel*FWR_Fiberlink.*OverlapMod.*dBLoss;
            obj.FiberLink(:,4) = OverlapMod*Fiber_Fresnel;
            
            
        end
        function[obj] = CombinedFP(obj)
            global curve_dir
            filename = 'detZyla42.mat';
            file = strcat(curve_dir,filename);
            load(file)
            obj.CombinedFPCh(:,1)= det_Zyla42(:,1);
            obj.CombinedFPCh(:,2)= det_Zyla42(:,2); 
            
        end
        
    end
    methods(Static)
        function[x_cut,y_cut] = select_bandpass(x,y,bounds)
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
        
    end
end

%Multiply Curves
function[x,y]= multiply_curves(x1,y1,x2,y2)

if x1(1)< x2(1)
    x = x1;
    Int = y1;
    wav2 = x2;
    Int2 = y2;
else
    wav2 = x1;
    Int2 = y1;
    x = x2;
    Int = y2;
end

[out1,~] = find(x == wav2(1),1);
x = x(out1:length(x));
Int = Int(out1:length(Int));

if x(end)> wav2(end)
else
    wav2_new = x;
    Int2_new = Int;
    
    wav_new = wav2;
    Int_new = Int2;
    x=wav_new;
    Int=Int_new;
    
    wav2=wav2_new;
    Int2=Int2_new;
end

[out1,~] = find(x == wav2(end),1,'last');

x = x(1:out1);
Int = Int(1:out1);

y = (Int).*(Int2);

end
%Dichroic curve maker
function [x,y] = Science_splitter(type)

%Science beamspliiter
if strcmp(type,'fiber')==1;
    t_width = 20;
    edges = [400 970 1310 2000];
    transmission = [5 95 5];
    reflect =0;
elseif strcmp(type,'guide')==1;
    t_width = 20;
    edges = [400 970 1300 2000];
    transmission = [5 95 5];
    reflect =1;
elseif strcmp(type,'LBTI')==1;
    t_width = 20;
    edges = [400 920 1300 2000];
    transmission = [5 95 95];
    reflect =0;
elseif strcmp(type,'AO')==1;
    t_width = 20;
    edges = [400 920 1300 2000];
    transmission = [5 95 95];
    reflect =1;
else
    disp('selected BS does not exist')
end

step = 1; % nm step size
pad{1}= [0 -t_width];
pad{2}= [0 0];
pad{3} = [t_width 0];
wave = [];
tave=[];
for ii = 1:3
    wavelength{ii} = (edges(ii)+pad{ii}(1):step:edges(ii+1)+pad{ii}(end));
    t{ii} = transmission(ii)*ones(size(wavelength{ii}));
    wave = [wave wavelength{ii}];
    tave = [tave t{ii}];
end

%% transition regions
transistion_wave{1} = (edges(2)-t_width+step:step:edges(2)-step);
transistion_wave{2} = (edges(3)+step:step:edges(3)+t_width-step);

%Interpolating with a spline shape
transition_wave = [transistion_wave{1} transistion_wave{2}];
transition_tave = interp1(wave,tave,transition_wave,'spline');

%% Combining specified + transition
[full_wave,I] = sort([wave transition_wave]);
full_tave = [tave transition_tave];
full_tave = full_tave(I);

% % check curve sensibility
% figure()
% plot(wave,tave,'.','markersize',7)
% figure(gcf)
% hold on
% plot(transition_wave,transition_tave,'.','markersize',7)
% plot(full_wave,full_tave+100,'.k','markersize',7)

full_tave=full_tave/100;

if reflect ==1;
    full_tave = 1-full_tave;
else
end
y = permute(full_tave,[2,1]);
x = permute(full_wave,[2,1]);





end
%Overlap Integral (Not used currently 4/7/17)
function [f,b,c,d,e] = Overlap_integral (e2_psf,e2_fiber,mu_PSF,mu_Fiber,dtheta)

% The coupling efficiency, q, is given in terms of the overlap integral of
% the field distributions of the incident beam and the fiber mode in
% the reference plane: integral(|Ef*Eb|)^2/ integral(|Ef|^2)integral(|Ef|^2)
% Source : Coupling of a semiconductor laser to a single-mode fiber
% Author : Nu Yu, February 1987

%all units are in microns!!

if length(e2_fiber) == 1;
    e2_fiber(1,2)=e2_fiber(1,1);
end

w1 = e2_psf(1,1)/2; %MFR (0.5/e^2) value of PSF
w2 = e2_psf(1,2)/2; %MFR (0.5/e^2) value of PSF

w3 = e2_fiber(1,1)/2; %MFR (0.5/e^2) value of fiber
w4 = e2_fiber(1,2)/2; %MFR (0.5/e^2) value of fiber

mu1 = mu_PSF(1,1); %offset of w1
mu2 = mu_PSF(1,2); %offset of w2
mu3 = mu_Fiber(1,1); %offset of w3
mu4 = mu_Fiber(1,2); %offset of w4

%initialize output arrays
b = [];
c = [];
d = [];
e = [];
%numerical method%

fun1 = @(x,y)(exp( -((x-mu1)/w1).^2).*exp( -((y-mu2)/w2).^2)).^2;  %generate 2D Gaussian function for PSF squared
fun2 = @(x,y)(exp( -((x-mu3)/w3).^2).*exp( -((y-mu4)/w4).^2)).^2;  %generate 2D Gaussian function for Fiber squared
fun3 = @(x,y) exp( -((x-mu1)/w1).^2).*exp( -((y-mu2)/w2).^2).*exp( -((x-mu3)/w3).^2).*exp( -((y-mu4)/w4).^2); %generate 2D Gaussian function for Fiber*PSF

%visualization of created functions
grid = 10*max(w1,w3); % create a dimension 10x larger than the waist
[X,Y]=meshgrid(-grid:grid,-grid:grid); % create a grid
PSF = fun1(X,Y); %make PSF
Fiber = fun2(X,Y);%make Fiber MFD
Total = PSF+Fiber;%combine on one grid
% imagesc(Total)

q1 = integral2(fun3,-inf,inf,-inf,inf); % integrate Fiber*PSF
q2 = integral2(fun2,-inf,inf,-inf,inf); % integrate Fiber
q3 = integral2(fun1,-inf,inf,-inf,inf); % integrate PSF
q4 = q1*q1; %square integral of Fiber*PSF
q5 = q4/(q2*q3);% normalize coupling

eta= (4*(w1^2*w3^2)/(w1^2+w3^2)^2)*exp(-2*((sqrt((mu1-mu3)^2+(mu2-mu4)^2))^2)/(w1^2+w3^2)); % check numerical method with simplified formula for photonics RP on fiber joints
q6 =(4*exp(-((2*mu1^2)/(w1^2 + w3^2)) + (4*mu1*mu3)/(w1^2 + w3^2) - (2*mu3^2)/(w1^2 + w3^2) - (2*(mu2 - mu4)^2)/(w2^2 + w4^2))...
    *sqrt(1/w1^2)*sqrt(1/w2^2)*sqrt(1/w3^2)*sqrt(1/w4^2))...
    /((1/w1^2 + 1/w3^2)*(1/w2^2 + 1/w4^2)); % derive the above formula using mathematica script (SMF_overlap_integral located in this folder) - unsimplified version

% For a circular PSF, with sigma_x = sigma_y and mu value, (b,c,d) will all give correct values for
% the overlap integral. % If elongated PSFs are used, then mehtods b and d
% must be used becasue c does not include 2 width parameters. For all test
% purposed so far c is not required but is left as a sanity check.


%angular tolerance
NA1 = 2*0.98/(pi*w3);
q7 = exp(-(2*dtheta/NA1)^2);

b = [b q5]; %q5 can be modified to work with Zemax PSF by performing numerical integral
c = [c eta];
d = [d q6];
e = [e q7];
f = b.*e;
end
%Wagner Static Losses
function [Strehl_static] = Wagner()
%--------------------------------------------------------------------------
% Gaussian Coupling Efficiency                                      General
% Description: Calculate Fiber coupling as found in Wagner and Tomlinson.
% Input:
% Variables:
% Sigma = Rms wavefront error
% Lambda = wavelength
% Omega = Wpp = Peak to Peak wavefront errror
% TR = Throughput losses due to random wavefront errors
% TR = exp(-(2*pi*sigma/lambda)^2)*T;
% Example: % Generic Throughput for random wave front perturbations
%--------------------------------------------------------------------------
T = 1;
lambda = 1;
% wpp = 0.633/10;%peak to valley (aka peak-peak)
wpp_gold = 0.633/20;%peak to valley (aka peak-peak)
wpp_die = 0.633/10;

rms_gold = wpp_gold/4;
num_optics_gold = 14;

rms_die = wpp_die/4;
num_optics_die = 0;

rms_total = sqrt(num_optics_gold*rms_gold^2+num_optics_die*rms_die^2);
TR = exp(-(2*pi*rms_total)^2)*T;
Strehl_static = TR;
end



