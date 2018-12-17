classdef Spectrograph < Instrument
    properties
        polarization
        grating
        maxR
        pixSamp
    end
    
    methods
        
        function [obj] = Spectrograph(type,polarization)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            %Need to define an optical model structure and rules.
            %Primary grating needs to be at the end
            
            % Any code not using output argument (obj)
            if nargin == 0
                
                type = 'Richardson';
                polarization = [0,0,0]; % degree of  polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)
                
            elseif nargin == 1
                
                type = 'Silicon';
                polarization = [0,0,0]; % degree of  polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)
                
            end
            
            
            if strcmp(type,'Old') == 1
                
                opticalModel{1} = struct('name','OffnerM1','type','OAP','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{2} = struct('name','OffnerM2','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{3} = struct('name','OffnerM3','type','OAP','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{4} = struct('name','OffnerM4','type','fold','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{5} = struct('name','OAP1','type','OAP','coatingName','FathomGold','number',2,'angle','10','efficiency',[],'polarization',0);
                opticalModel{6} = struct('name','Spectrum Mirror','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{7} = struct('name','OAP2','type','OAP','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{8} = struct('name','Cross Disperser','type','grating','coatingName','Xdisp','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{9} = struct('name','TMA1','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{10} = struct('name','TMA2','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{11} = struct('name','TMA3','type','spheric','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{12} = struct('name','H4RG','type','detector','coatingName','detH4RG','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{13} = struct('name','R6Grating','type','Grating','coatingName','R6Grating','number',1,'angle','81','efficiency',[],'polarization',0);
                bandPass = [970,1270];
                
            elseif strcmp(type,'Richardson') == 1
                
                opticalModel{1} = struct('name','Col M1','type','conic','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{2} = struct('name','Col M2','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{3} = struct('name','Col M3','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{4} = struct('name','Fold M1','type','fold','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{5} = struct('name','Cross Disperser','type','grating','coatingName','Xdisp','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{6} = struct('name','Cam M1','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{7} = struct('name','Cam M2','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{8} = struct('name','Cam M3','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{9} = struct('name','H4RG','type','detector','coatingName','detH4RG','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{10} = struct('name','R6Grating','type','Grating','coatingName','R6Grating','number',1,'angle','81','efficiency',[],'polarization',0);
                bandPass = [800,1600];
                
                
            elseif strcmp(type,'Silicon') == 1
                
                opticalModel{1} = struct('name','Col M1','type','conic','coatingName','FathomGold','number',1,'angle','10','efficiency',[],'polarization',0);
                opticalModel{2} = struct('name','Col M2','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{3} = struct('name','Col M3','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{4} = struct('name','Fold M1','type','fold','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{5} = struct('name','Cross Disperser','type','grating','coatingName','Xdisp','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{6} = struct('name','Cam M1','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{7} = struct('name','Cam M2','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{8} = struct('name','Cam M3','type','conic','coatingName','FathomGold','number',1,'angle','25','efficiency',[],'polarization',0);
                opticalModel{9} = struct('name','H4RG','type','detector','coatingName','detH4RG','number',1,'angle',[],'efficiency',[],'polarization',0);
                opticalModel{10} = struct('name','R6Grating','type','Grating','coatingName','R6Silicon','number',1,'angle','81','efficiency',[],'polarization',0);
                bandPass = [800,1600];
            end
            
            
            %             polarization = [0.5,0.1,0]; % degree of  polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)
            
            current_path = pwd;
            if strcmp(current_path(2:8),'Volumes')==1
                
                curveDirectory = '/Volumes/Software/Simulator/RefFiles/Curves/Spectrograph/';
                load('/Volumes/Software/Simulator/polycoeffs2.mat')
            elseif strcmp(current_path(2:4),'afs')==1
                
                curveDirectory = '/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/RefFiles/Curves/Spectrograph/';
                load('/afs/crc.nd.edu/group/Exoplanets/ebechter/NewSim/Simulator/polycoeffs2.mat')               
            else
                
                curveDirectory = [current_path(1:2) '\Simulator\RefFiles\Curves\Spectrograph\'];
                load([current_path(1:2) '\Simulator\polycoeffs2.mat']);
            end
            
            
            
            maxR = 275e3;
            pixSamp = 3;
            
            
            obj.bandPass = bandPass;
            obj.polarization = polarization;
            obj.opticalModel = opticalModel;
            obj.maxR = maxR;
            obj.pixSamp = pixSamp;
            obj.name = 'Spectrograph';
            [obj] = loadOpticalModelCurves(obj,curveDirectory);
            [obj] = trimThroughput(obj);
            [obj] = R6Grating(obj,curveDirectory);
            [obj] = Include_Grating(obj);
            
            
            detBound = [-20.48:1:+20.48]; %mm centered (slightly over the detector)
            
            for ii = 1:size(wave_coeff,1)
                det_wave = polyval(wave_coeff(ii,:,1),detBound);
                det_eff(:,ii) = interp1(obj.finalThroughput{1,2}(:,ii),obj.finalThroughput{1,1}(:,ii),det_wave*1000);
            end
            
            obj.intTrans = mean(mean(det_eff));
        end
        
        function[obj] = R6Grating(obj,curveDirectory)
            %% ----------Load Curves----------%
            gratingfile = [curveDirectory,'R6Grating81.mat'];
            load(gratingfile)
            
            %% ----------Assign Object Properties----------%
            if isempty (obj.polarization == 1) || obj.polarization(1,3) == 0 || obj.opticalModel{ii}.polarization == 0
                % if the user does not set the polarization, sets it but
                % flags as 0 or turns off the grating polarization
                % specifically
                % obj.grating= GratingEff_new;
                % !!!!!!!!!at PDR we used this but now we should use will use
                % unpolarized grating set below.
                pfrac = 0;
                sfrac = 1-pfrac;
                dop =  0; % this makes it unpolarized
                
            else
                % if polarization state is set, is flagged as 1 and the
                % grating polarization is turned on
                pfrac = obj.polarization(1,2);
                sfrac = 1-pfrac;
                dop =  obj.polarization(1,1);
            end
            
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
            obj.grating(:,1)=GratingEff_new(:,1);
            obj.grating(:,2:40)=polarized;
            
        end
        
        function [obj] = Include_Grating(obj)
            GratingEff = obj.grating;
            for ii = 2:size(GratingEff-3,2) % trying to fix for 36 orders
                y1 = obj.finalThroughput(:,2);
                x1 = obj.finalThroughput(:,1);
                y2 = GratingEff(:,ii);
                x2 = GratingEff(:,1);
                yq = interp1(x1,y1,x2,'linear','extrap');
                [wav_order,Tput_order]= Instrument.multiply_curves(x2,y2,x2,yq);
                clear xq vq
                orders{1}(:,ii-1)=Tput_order;
                orders{2}(:,ii-1)=wav_order;
                orders{3}(ii-1)=156-ii-1;
                
            end
            %----------Assign Object properties----------%
            obj.finalThroughput = orders;
        end
        
        function [] = spectrographPlot(obj)
            handle =[];
            
            %% Custom color lists, yo
            d = get(groot,'DefaultAxesColorOrder');
            for ii = 1:7
                colors{ii}=d(ii,:);
            end
            colors{8}= [0.175 0.175 0.175];
            colors{9}= colors{2};
            colors{10}= colors{3};
            colors{11} =[0 0.3 0];
            colors{12} =colors{4};
            colors{13} =colors{5};
            clear d
            
            labels = [];
            for ii = 1:size(obj.opticalModel,2)
                temp{1} = obj.opticalModel{ii}.name;
                labels = [labels;temp];
            end
            
            figure
            hold on
            for ii = 1:length(obj.progress)
                h{ii} = plot(obj.progress{ii}(:,1),obj.progress{ii}(:,2),'.','Markersize',8,'Color',colors{ii});
                handle = [handle h{ii}];
            end
            
            for jj = 1:size(obj.finalThroughput{1},2)
                h{ii+1}=plot(obj.finalThroughput{2}(:,jj),obj.finalThroughput{1}(:,jj),'.','Color',colors{11},'Markersize',8);
                if jj == 1
                    handle = [handle h{ii+1}];
                end
            end
            l=legend(handle,labels,'Location','best');
            plot_max = max(obj.progress{length(obj.progress)}(:,1));
            plot_min = min(obj.progress{length(obj.progress)}(:,1));
            ylim([0 1])
            xlim([900 1350])
            ylabel('Throughput')
            xlabel('\lambda nm')
            l.FontSize = 10;
            l.Box = 'off';
            box on
            ax = gca;
            ax.LineWidth = 1.5;
        end
    end
end