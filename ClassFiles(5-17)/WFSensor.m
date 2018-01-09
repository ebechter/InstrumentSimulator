classdef WFSensor < Instrument
    properties
        EntWindowWFS
        WFS
        AO
    end
    methods
        function[obj] = WFSensor(AO,entwindow)
            obj.AO = AO;
            
            if nargin == 0
                % Provide values for superclass constructor
                % and initialize other inputs
                 AO = 1;
                 obj.Type = 'LBTI AO';
            else
                % need to put conditions for user specified inputs here
                if AO ==1
                    obj.Type = 'LBTI AO';
                elseif AO == 2
                    obj.Type = 'SOUL 1"';
                elseif AO == 3
                    obj.Type = 'SOUL 0.8"';
                elseif AO == 4
                    obj.Type = 'SOUL 0.6"';
                end
            end
            obj = LBTMirrors(obj);
            obj = EntranceWindowWFS(obj,entwindow);
            obj = WFSOptics(obj,AO);
            obj.PathName = 'AO'; %Choose Path Name 'Fiber','Image','Quad','AO','Spectrograph','WFC'
            obj = Optical_Path(obj); %Set Optical Path parameters
            obj = Path_Multiply(obj); %Mulitply curves
            obj = Trim_Throughput(obj);%Trim to Bandpass (1,:)
            obj = Integrated_Transmission(obj); %Calc Integrated Throughput
        end
        function[obj] = EntranceWindowWFS(obj,entwindow)
            global curve_dir
            
            if strcmp(entwindow,'ECI') ==1 
            %% ----------Load Curves----------%
            filename = 'EntWRefl.mat';
            file = strcat(curve_dir,filename);
            load(file)
            LBTI_bs_wave = EntWRefl(:,1);
            LBTI_bs_trans = EntWRefl(:,2);
            
            elseif strcmp(entwindow,'sim') ==1 
            %% ----------Generate Curves----------%
            type = 'AO';
            [LBTI_bs_wave,LBTI_bs_trans]=Science_splitter(type);
            
            elseif strcmp(entwindow,'existing') ==1
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
            
            %----------Load Curves----------%
            filename = 'detOCAM.mat';
            file = strcat(curve_dir,filename);
            load(file)
            
            filename = 'detCCD39.mat';
            file = strcat(curve_dir,filename);
            load(file)
            %----------Multiply Curves----------%
            FLAO = Optics.* det_CCD39(:,2);
            SOUL = Optics.* det_OCAM(:,2);
            %----------Assign Object properties----------%
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
    end
end

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





end %Dichroic curve maker
