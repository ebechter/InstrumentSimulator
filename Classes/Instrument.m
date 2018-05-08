classdef Instrument
    
    properties
        opticalModel
        name
        bandPass
        progress
        finalThroughput
        intTrans
    end
    
    methods
        
        function [obj] = loadOpticalModelCurves(obj,curveDirectory)
            
            if strcmp(class(obj),'Spectrograph') == 1
                n = size(obj.opticalModel,2)-1; %number of optics minus primary grating
            else
                n = size(obj.opticalModel,2); %number of optics 
            end
            
            for ii =1:n %loop over cellarry size (i.e. surfaces)

                number = obj.opticalModel{ii}.number;
                angle = obj.opticalModel{ii}.angle;
                coatingName = obj.opticalModel{ii}.coatingName;
                filename = strcat(curveDirectory,coatingName,angle);
                
           
                
                if isa(obj,'Imager') || isa(obj,'AO') || ...
                        isempty (obj.opticalModel{ii}.polarization == 1) || isempty (obj.polarization == 1) || obj.polarization(1,3) == 0
                    %if the user does not want polarization to be used or
                    %is not set
                    pol = 0; % set pol = 0 to use average curve
                    obj.opticalModel{ii}.efficiency = Instrument.loadCurves(filename,number,coatingName);
                else
                    %if the user does want polarization to be used
                    pol = obj.opticalModel{ii}.polarization;
                    
                    pfrac = obj.polarization(1,2);
                    dop =  obj.polarization(1,1);
                    
                    obj.opticalModel{ii}.efficiency = Instrument.loadPolCurves(filename,number,coatingName,pfrac,dop);
                end
                
                
                if ii == 1
                    pgress{1} = obj.opticalModel{1}.efficiency;
                else
                    [pgress{ii}(:,1),pgress{ii}(:,2)]= Instrument.multiply_curves(pgress{ii-1}(:,1),pgress{ii-1}(:,2),obj.opticalModel{ii}.efficiency(:,1),obj.opticalModel{ii}.efficiency(:,2));
                end
            end
            
            obj.progress = pgress;
        end
        function [obj] = trimThroughput(obj)
            for ii = 1:length(obj.progress)
                [xtrim,ytrim]=Instrument.trimToBand(obj.progress{ii}(:,1),obj.progress{ii}(:,2),obj.bandPass(1,:));
                obj.progress{ii} = [xtrim,ytrim];
            end
            obj.finalThroughput(:,1) = xtrim;
            obj.finalThroughput(:,2) = ytrim;
        end  
        function [] = progressPlot(obj)
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

            l=legend(handle,labels,'Location','best');
            plot_max = max(obj.progress{length(obj.progress)}(:,1));
            plot_min = min(obj.progress{length(obj.progress)}(:,1));
            ylim([0 1])
            xlim([plot_min plot_max])
            ylabel('Throughput')
            xlabel('\lambda nm')
            l.FontSize = 10;
            l.Box = 'off';
            box on
            
        end        
        function[obj] = intThroughput(obj)
            total = (max(obj.bandPass(1,:))-min(obj.bandPass(1,:)));
            obj.intTrans= trapz(obj.finalThroughput(:,1),obj.finalThroughput(:,2))/total;
        end
    end
    
    methods(Static)
        function [throughput] = loadCurves(filename,number,coatingName)
            
            temp = load([filename,'.mat']); % convert .mat to normal array from default struct format
            temp = struct2cell(temp);
            surface = temp{1};
            
            if strcmp(coatingName,'FathomGold')==1
                % WL(nm),10deg r-pol,10deg s-pol,10deg p-pol,25deg r-pol,25deg s-pol,25deg p-pol
                Spol = (surface(:,3)/100);
                Ppol = (surface(:,4)/100);
                Rpol = (Spol+Ppol)/2;
            else
                Rpol = surface(:,2);
            end
            throughput(:,1) = surface(:,1);
            throughput(:,2) = Rpol.^number;
        end
        function [throughput] = loadPolCurves(filename,number,coatingName,pfrac,dop)
            
            temp = load([filename,'.mat']); % convert .mat to normal array from default struct format
            temp = struct2cell(temp);
            surface = temp{1};
            
            if strcmp(coatingName,'FathomGold')==1
                Spol = (surface(:,3)/100);
                Ppol = (surface(:,4)/100);
                sfrac = 1-pfrac;
                unpolarized = 0.5*(Spol+Ppol);
                
                Rpol = dop*(Spol*sfrac+Ppol*pfrac)+unpolarized*(1-dop);                
                %Rpol = (Spol+Ppol)/2;
            else
                Rpol = surface(:,2);
            end
            throughput(:,1) = surface(:,1);
            throughput(:,2) = Rpol.^number;
        end
        
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
        function[x_cut,y_cut] = trimToBand(x,y, bounds)
            rng_l =max(bounds(:,1));
            rng_u =min(bounds(:,2));
            wavelength_limits = x<=rng_u & x>=rng_l; % Instrument band
            y_cut = y(wavelength_limits);% cut the flux value at the band
            x_cut = x(wavelength_limits);
        end
    end
end





