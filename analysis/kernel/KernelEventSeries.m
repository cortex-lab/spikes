

classdef KernelEventSeries < handle
    
    properties
        
        eventTimes
        eventValues
        
        eventName
        window
        
        predictors = [];
        
        upSampFactor = 1; % needs to match KernelModel's
    end
    
    methods
        
        function obj = KernelEventSeries(eventName, eventTimes, eventValues, window)
            obj.eventName = eventName;
            obj.eventTimes = eventTimes;
            obj.eventValues = eventValues;
            obj.window = window;
            
        end
        
        function t = predictorTimes(obj, timeBinSize)
            nSampExpected = obj.upSampFactor*ceil(diff(obj.window)/timeBinSize);
            t = obj.window(1)+(0:nSampExpected-1)*timeBinSize/obj.upSampFactor;
        end
        
        function addPredictor(obj, pr, timeBinSize)
            % the predictor is a vector which must cover the obj.window
            % with time resolution of 10x the model's time bin size
            
            nSampExpected = numel(obj.predictorTimes(timeBinSize));
            if numel(pr)~=nSampExpected
                error('Wrong number of samples in the predictor');                
            end
            
            if isempty(obj.predictors)
                obj.predictors = pr(:);
            else
                obj.predictors(:,end+1) = pr;
            end
        end
        
    end
    
end
