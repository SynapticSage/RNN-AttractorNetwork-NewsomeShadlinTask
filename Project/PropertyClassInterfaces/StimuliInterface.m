classdef StimuliInterface % Abstract Class
    % Template class for organizing the parameters and computation of
    % stimuli that will be applied in the simulation.
    
    % INPUT
    properties (Access = public)
        I_base = 1;     % Default scale of the stimulus .. 
        neurIdentities; % Universal identity bector that labels how many neurons of what types, and in which locations in the vector
    end
    
    % Methods that determine stimulus properties
    methods (Abstract)
        % Abstract method that handles that creates stimuli and weights
        % across the neural population, given correct parameters. This
        % method can be implemented for any stimulus regime.
        generateStimuli(this); 
    end
    
    % OUTPUT
    properties (SetAccess = protected)
        I_W; % Weights to the stimulus, uniform vector if all units receive same input
        Iapp;
    end
    
end