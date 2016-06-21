classdef Neuron < handle
    
    % Object Properties
    properties
    end
    
    % Biophysical Properties
    properties (SetObservable = true, Access = public)
        
        rmax;       % max firing rate
        rthresh;    % threshold
        rsigma;     % how currents chage to cell
        
        Synapses;   % set of objects containing information about synapses
    end
    
    methods (Abstract)
        
    end
    
end