classdef Spectrograph < Instrument
    properties
        polarization
        grating
        maxR
        pixSamp
    end
    
    methods
        function [obj] = Spectrograph(bandPass, polarization, opticalModel,maxR,pixSamp)
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                %Need to define an optical model structure and rules.
                %Primary grating needs to be at the end
                opticalModel{1} = struct('name','OffnerM1','type','fold','coatingName','FathomGold','number',2,'angle','10','efficiency',[]);
                opticalModel{2} = struct('name','OffnerM2','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[]);
                opticalModel{3} = struct('name','FoldM1','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[]);
                opticalModel{4} = struct('name','M2','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[]);
                opticalModel{5} = struct('name','M3','type','fold','coatingName','FathomGold','number',1,'angle','10','efficiency',[]);
                opticalModel{6} = struct('name','R6Grating','type','Grating','coatingName','R6Grating','number',1,'angle','88','efficiency',[]);
                bandPass = [970,1270];
                polarization = [0.5,0.1,0]; % degree of polarization, P-fraction, flag (1 has pol effects, 0 reverts to original)
                curveDirectory = [pwd '/RefFiles/Curves/Spectrograph/'];                
                maxR = 275e3;
                pixSamp = 3;
                
            else
                
            end
            
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
        end
        
        function[obj] = R6Grating(obj,curveDirectory)
            %% ----------Load Curves----------%
           gratingfile = [curveDirectory,'R6Grating81.mat'];
           load(gratingfile) 
           
            %% ----------Assign Object Properties----------%
            if isempty (obj.polarization == 1) || obj.polarization(1,3) == 0
            obj.grating= GratingEff_new;
            else
               pfrac = obj.polarization(1,2);
               sfrac = 1-pfrac;
               dop =  obj.polarization(1,1);
               
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