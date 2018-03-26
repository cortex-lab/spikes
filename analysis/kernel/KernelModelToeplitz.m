
classdef KernelModelToeplitz < KernelModel
    
    properties
        
    end
    
    methods
        
        function obj = KernelModelToeplitz(tN, lambda)
            obj@KernelModel(tN, lambda);
        end
        
        function generatePredictor(obj)
            % generate the predictor matrix
            
            % determine number of predictors
            
            % for each event
            for es = 1:numel(obj.eventSeries)
                
                e = obj.eventSeries{e};
                
                
            
        end
        
    end
    
end