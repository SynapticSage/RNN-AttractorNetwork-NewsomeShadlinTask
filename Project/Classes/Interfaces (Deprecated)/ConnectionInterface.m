classdef ConnectionInterface % Abstract Class
    
    % General class the handles the act of assigning connection between
    % neurons and neuronal groups. The type of connections implemented
    % depend on the implementation of this class.
    
    % INPUTS
    properties (Access = public)
        
    end
    
    % Methods that control assignment of the W given the rules of the class
    methods (Abstract)
        generateConnections(this);    % abstract method for controlling the specification of weights, dependent on neural subtypes
        returnOutputs(this);
    end
    
    % OUTPUTS
    properties (Access = public)
        
    end
end