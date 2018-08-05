

classdef KernelModel < handle
    
    properties
        name
        timeBinSize
        eventSeries % objects that 
        lambda = 0; % regularization term
        A = []; % predictors   
        X = [];
        tN = []; 
        nPred = 0;
        predictorInds = [];
        upSampFactor = 1;
    end
    
    methods
        function obj = KernelModel(name, tN, lambda)
            obj.timeBinSize = mean(diff(tN)); % assuming evenly sampled in time
            obj.lambda = lambda;
            obj.tN = tN;
            obj.name = name;
        end
                
        function addEventSeries(obj, es)
            % add a new set of events to regress
            obj.eventSeries{end+1} = es;            
        end
        
        function setLambda(obj, lambda)
            % change the normalization value
            obj.lambda = lambda;            
            if isempty(lambda) && ~isempty(obj.A)
                % remove normalization
                obj.A = obj.A(1:numel(obj.tN),:);
            elseif ~isempty(lambda)
                % set the normalization                
                obj.A(numel(obj.tN)+1:numel(obj.tN)+obj.nPred, 1:end-1) = ...
                    diag(lambda*ones(1,obj.nPred));
            end
            
        end
        
        function generatePredictor(obj)
            
            % determine number of predictors
            np = 0;
            predictorInds = [];
            for es = 1:numel(obj.eventSeries)
                e = obj.eventSeries{es};
                np = np+size(e.predictors,2);
                predictorInds(end+1:end+size(e.predictors,2)) = es;
            end
            obj.predictorInds = predictorInds;
            obj.nPred = np;
            A = zeros(numel(obj.tN), np);            
            
            % for each event
            pInd = 1;
            for es = 1:numel(obj.eventSeries)
                
                e = obj.eventSeries{es};
                thisEv = zeros(1,numel(obj.tN)*obj.upSampFactor);
                evInds = round((e.eventTimes-obj.tN(1)+e.window(1))/obj.timeBinSize*obj.upSampFactor);
                thisEv(evInds) = e.eventValues;
                thisT = obj.tN(1)+(0:numel(thisEv)-1)*obj.timeBinSize/obj.upSampFactor;
                
                for p = 1:size(e.predictors,2)
                    aa = conv(thisEv, e.predictors(:,p), 'full');                    
                    A(:,pInd) = interp1(thisT, aa(1:numel(thisEv)), obj.tN);
                    pInd = pInd+1;
                end
            end
                    
            A(:,end+1) = 1; % intercept term
            
            obj.A = A;
            
            obj.setLambda(obj.lambda);
        end
        
        function [k, t] = kernelFor(obj, eventName)
            % return the kernel for an event: multiply the fit by the basis
            % functions
            eventInd = find(cellfun(@(x)strcmp(eventName, x.eventName), obj.eventSeries));
            if isempty(eventInd)
                warning('could not find event %s', eventName);
            else
                e = obj.eventSeries{eventInd};
                d = e.predictors;
                t = e.predictorTimes(obj.timeBinSize);                
                k = d*obj.X(obj.predictorInds==eventInd,:);
            end
        end
        
        function pred = predictData(obj)
            if isempty(obj.X)
                error('add the results of a fit (X) first');
            end
                        
            pred = obj.A(1:numel(obj.tN),:)*obj.X;
        
        end
        
    end
    
end
        