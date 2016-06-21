classdef Simulation < handle
    % Class Name:   Solver
    % Pupose:       Handles the algorithmic portion of the differential
    %               equations setup by the neurons and networking.
    %               Chooising an object-oriented approach because I'm
    %               trying to create a cleaner, easier to maintain, more
    %               plug-and-play style simulation.
    
    % Simulation temporal precision
    properties
    end
    
    % Differential equation system variables
    properties
        I;          % Record of all current input states
        r;          % Record of all current firing rates
        D;          % Record of all depression states
        S;          % Record of all synaptic transmission states
    end
    
    % Tab keepers on objects that enter simulation
    properties
    end
    
    methods
        function this = Simulation(~)
            
        end
    end
    
end