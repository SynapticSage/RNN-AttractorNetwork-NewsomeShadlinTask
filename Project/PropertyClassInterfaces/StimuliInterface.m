classdef StimuliInterface % Abstract Class
    % Template class for organizing the parameters and computation of
    % stimuli that will be applied in the simulation.
    
    % INPUT
    properties (Access = public)
        % General stimulus property
        I_base = 1;     % Default scale of the stimulus .. 
        
        % General things that the stimulus has to have acces to to
        % calculate inputs wrt neurons and model time
        neurIdentities; % Universal identity bector that labels how many neurons of what types, and in which locations in the vector
        dt;
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
        % General applied current due to the stimulus to the n units in the
        % simulation.
         Iapp;
    end
    
end