classdef ConnectionInterface % Abstract Class
    
    % General class the handles the act of assigning connection between
    % neurons and neuronal groups. The type of connections implemented
    % depend on the implementation of this class.
    
    % INPUTS
    properties (Access = public)
        neurIdentities; % Stores information about the neural units/groups, how they present across the vector -- it's how we label which neurons are which types, and that determines by the rules in the assignment method how to setup the rules.
    end
    
    % Methods that control assignment of the W given the rules of the class
    methods (Abstract)
        generateConnections(this);    % abstract method for controlling the specification of weights, dependent on neural subtypes
        returnOutputs(this);
    end
    
    % OUTPUTS
    properties (Access = public)
        W; % Stores the resulting connection strenghts
    end
end