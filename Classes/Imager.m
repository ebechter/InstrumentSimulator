classdef Imager < Instrument
    properties
        psf
        pixelPitch
        detectorDimensions
        plateScale
        apDiameter
        blockFrac
    end
    
    methods
        function [obj] = Imager(type,bandPass,opticalModel,pixelPitch,detectorDimensions,psf)
            %% Pre Initialization %%
            
            blockFrac = [];
            apDiameter = [];
            % Any code not using output argument (obj)
            if nargin == 0 || strcmp(type,'ANDOR') == 1
                %Need to define an optical model structure and rules.
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','L1','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{5} = struct('name','ADC1','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{6} = struct('name','ADC2','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{7} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{8} = struct('name','Shortpass','type','triplet prism','coatingName','ECISpassT','number',1,'angle','45','efficiency',[]);
                %opticalModel{9} = struct('name','AndorFilter','type','bandpass filter','coatingName','Imagefilter','number',1,'angle',[],'efficiency',[]);
                opticalModel{9}= struct('name','L5','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{10}= struct('name','ANDOR','type','detector','coatingName','AndorZyla42','number',1,'angle',[],'efficiency',[]);
                
                plateScale = [];
                pixelPitch = [];
                detectorDimensions = [];
                bandPass = [900,970];
                
                
            elseif strcmp(type,'WFC') == 1
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','WFCL1','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{5} = struct('name','WFCL2','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{6} = struct('name','Basler','type','detector','coatingName','Basler','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [600,1200];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale =[];
            
             elseif strcmp(type,'FIBER') == 1
                
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','L1','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{5} = struct('name','ADC1','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{6} = struct('name','ADC2','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{7} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{8} = struct('name','Shortpass','type','Dichroic','coatingName','ECISpassR','number',1,'angle','45','efficiency',[]);
                opticalModel{9} = struct('name','Longpass','type','Dichroic','coatingName','ECILpassR','number',1,'angle','45','efficiency',[]);
                opticalModel{10}= struct('name','L2','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{11}= struct('name','L3','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                
               
                bandPass = [970,1300];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                
            elseif strcmp(type,'QUADCELL') == 1
                
                opticalModel{1} = struct('name','M1','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{2} = struct('name','M2','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','M3','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','L1','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{5} = struct('name','ADC1','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{6} = struct('name','ADC2','type','triplet prism','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{7} = struct('name','FSM','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{8} = struct('name','Shortpass','type','Dichroic','coatingName','ECISpassR','number',1,'angle','45','efficiency',[]);
                opticalModel{9} = struct('name','Longpass','type','Dichroic','coatingName','ECILpassT','number',1,'angle','45','efficiency',[]);
                opticalModel{10}= struct('name','L5','type','triplet lens','coatingName','NIR','number',1,'angle',[],'efficiency',[]);
                opticalModel{11}= struct('name','QuadCell','type','detector','coatingName','NIRQuadCell','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [1000,1500];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                           
            elseif strcmp(type,'LBT') == 1
                opticalModel{1} = struct('name','Primary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
                opticalModel{2} = struct('name','Secondary','type','mirror','coatingName','LBT_T_Al','number',1,'angle',[],'efficiency',[]);
                opticalModel{3} = struct('name','Tertiary','type','mirror','coatingName','LBT_PS_Al','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [400,2000];
                pixelPitch = [];
                detectorDimensions = [];
                plateScale = [];
                blockFrac = 0.108; % value listed as 0.889
                apDiameter = 8.22; % meters (IR effective collecting area, supposedly)
                
                
            elseif strcmp(type,'LBTI') == 1
                opticalModel{1} = struct('name','Entrance Window','type','dichroic','coatingName','ECIEntWT','number',1,'angle',[],'efficiency',[]);
                opticalModel{2} = struct('name','Ellipse Mirror','type','elliptical mirror','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{3} = struct('name','Pupil mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{4} = struct('name','Roof mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{5} = struct('name','Switch Box mirrror','type','fold','coatingName','ProtectedGold','number',1,'angle','45','efficiency',[]);
                opticalModel{6} = struct('name','Exit window','type','mirror','coatingName','SilicaAR','number',1,'angle',[],'efficiency',[]);
                
                bandPass = [400,2000];
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
