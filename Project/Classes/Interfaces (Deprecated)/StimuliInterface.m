classdef StimuliInterface % Abstract Class
    % Template class for organizing the parameters and computation of
    % stimuli that will be applied in the simulation.
    
    % INPUT
    properties (Access = public)
     
    end
    
    % Methods that determine stimulus properties
    methods (Abstract)
        % Abstract method that handles that creates stimuli and weights
        % across the neural population, given correct parameters. This
        % method can be implemented for any stimulus regime.
        generateStimuli(this); 
        returnOutputs(this);
    end
    
    % OUTPUT
    properties (SetAccess = protected)
        
    end
    
end