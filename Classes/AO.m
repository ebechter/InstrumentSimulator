classdef AO < Instrument
    properties
        psf
        pixelPitch
        detectorDimensions
        plateScale
        apDiameter
        blockFrac
    end
    
    methods
        function [obj] = AO(type,bandPass,opticalModel,pixelPitch,detectorDimensions,psf)
            %% Pre Initialization %%
            
            blockFrac = [];
            apDiameter = [];
            % Any code not using output argument (obj)
            
            opticalModel{1} = struct('name','Primary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
            opticalModel{2} = struct('name','Secondary','type','mirror','coatingName','LBT_T_Al','number',1,'angle',[],'efficiency',[]);
            opticalModel{3} = struct('name','Tertiary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
            
            if nargin == 0
                
                % use default FLAO settings
                
                opticalModel{4} = struct('name','EntWindow','type','dichroic','coatingName','LBTI_ZnSeR','number',1,'angle',[],'efficiency',[]);
                opticalModel{5} = struct('name','InternalOptics','type','single number','coatingName','WFS_Internal','number',1,'angle',[],'efficiency',[]); %0.9635^11
                opticalModel{6} = struct('name','CCD39','type','detector','coatingName','CCD39','number',1,'angle',[],'efficiency',[]);
                
                
                
            elseif strcmp(type(1:4),'FLAO') == 1
                
                if strcmp(type(5:end),'ilocater') == 1
                    opticalModel{4} = struct('name','EntWindow','type','dichroic','coatingName','ECIEntWR','number',1,'angle',[],'efficiency',[]);
                else
                    opticalModel{4} = struct('name','EntWindow','type','dichroic','coatingName','LBTI_ZnSeR','number',1,'angle',[],'efficiency',[]);
                end
                
                opticalModel{5} = struct('name','InternalOptics','type','single number','coatingName','WFS_Internal','number',1,'angle',[],'efficiency',[]); %0.9635^11
                opticalModel{6} = struct('name','CCD39','type','detector','coatingName','CCD39','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [600,1000];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                
            elseif strcmp(type(1:4),'SOUL') == 1
                
                if strcmp(type(5:end),'ilocater') == 1
                    opticalModel{4} = struct('name','Entrance Window','type','dichroic','coatingName','ECIEntWR','number',1,'angle',[],'efficiency',[]);
                else
                    opticalModel{4} = struct('name','Entrance Window','type','dichroic','coatingName','LBTI_ZnSeR','number',1,'angle',[],'efficiency',[]);
                end
                opticalModel{5} = struct('name','InternalOptics','type','single number','coatingName','WFS_Internal','number',1,'angle',[],'efficiency',[]); %0.9635^11
                opticalModel{6} = struct('name','OCAM2K','type','detector','coatingName','OCAM2K','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [600,900];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                
            end
            curveDirectory = [pwd '/RefFiles/Curves/Instrument/'];
            obj.bandPass = bandPass;
            obj.pixelPitch = [];
            obj.detectorDimensions = [];
            obj.opticalModel = opticalModel;
            obj.plateScale = plateScale;
            obj.blockFrac = blockFrac;
            obj.apDiameter = apDiameter;
            [obj] = loadOpticalModelCurves(obj,curveDirectory);
            [obj] = trimThroughput(obj);
            [obj] = intThroughput(obj);
        end
    end
    
    methods(Static)
        
    end
end
