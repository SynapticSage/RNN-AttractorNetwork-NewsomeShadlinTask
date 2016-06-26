classdef NeuronProperties
    
    properties
        % Vector marking how many excitatory and how many inhibitory
        % neurons as well as their (arbitrary) positions in a vector, which
        % subsequently determines how other vectors populate
        neurIdentities;
        
        % Maximum of the excitatory and inhbitory neurons
        rMax0E;
        rMax0I;
        
        % Specification for release probabilities
        p0E;
        p0I;
        
        % Specifications for synaptic time constants
        tausE;
        tausI;
        
        % Specification for depression time constants
        tauDbar;
        tauDvar=0;
        
        % Specification for membrane time constants
        tauMbar;
        tauMvar=0;
        
        % Defautl response to input stimuli
        Ith0E;
        Ith0I;
        
        % Input width
        IwidthE;
        IwidthI;
        
        % Isigma threshold
        IsigmaThE=0;
        IsigmaThI=0;
        
        % Synaptic release/binding probs/fractions
        sFrac=1;

        % whether use gpuArray in output
        useGPU=false;
       
    end
    
    properties (Access=protected)
        initFunction; % Stores the initialization function used to initialize output variables. This controls whether or not it creates gpu arrays
    end
    
    methods
        % ---------------------------------------------------------------
        function this = NeuronProperties(this) %#ok<INUSD>
            % Constructor method that primarily determine how the generator
            % should initialize it's data types, gpu-driven or not
            
            fprintf('Please set neuron properties ...\n');
            
        end
        % ---------------------------------------------------------------
        function this= generateOutputParams(this)
            % This function takes all of the parameters given at the outset
            % of the model related to connnections and generates the vector
            % formats needed by the model loop
            
            % Create an alias for this object
            t=this;
            
            % Create a rand stream for any random properties of the vectors
            % created per neuron
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');
            
            % Get indices of the exitatory and inhibitory cells
            exc         = t.neurIdentities == 0;
            inh         = t.neurIdentities == 1;
            
            % Assign the rMax vector
            t.rMax = zeros(size(t.neurIdentities));
            t.rMax(exc) = t.rMax0E;
            t.rMax(inh) = t.rMax0I;
            
            % Assign the synaptic time const
            t.tauS = zeros(size(t.neurIdentities));
            t.tauS(exc) = t.tausE;
            t.tauS(inh) = t.tausI;
            
            % Assign the depression time const...
            % creates a random constant for the excitatory cells, and a
            % non-random for inhibitory cells
            t.tauD = zeros(size(t.neurIdentities));
            t.tauD(exc) = t.tauDbar + t.tauDvar*rand(r,size(t.tauD(exc)));
            t.tauD(inh) = t.tauDbar;
            
            % Assign the membrane time constants
            t.tauM = zeros(size(t.neurIdentities));
            t.tauM(exc) = t.tauMbar + t.tauMvar*rand(r,size(t.tauM(exc)));
            t.tauM(inh) = t.tauMbar;
            
            % Specify the release propabilities per synapse
            t.p0 = zeros(size(t.neurIdentities));
            t.p0(exc) = t.p0E;
            t.p0(inh) = t.p0I;
            
            % Specify the threshold input and input width
            t.Iwidth    = zeros(size(t.neurIdentities));
            t.Ith       = zeros(size(t.neurIdentities));
            t.Iwidth(exc) = t.IwidthE;
            t.Iwidth(inh) = t.IwidthI;
            t.Ith(exc)  = t.Ith0E + rand(size(t.Ith(exc)))*t.IsigmaThE;
            t.Ith(inh)  = t.Ith0I + rand(size(t.Ith(inh)))*t.IsigmaThI;
            
            % Specify the sfrac
            if isempty(t.sFrac)
                error('Input sFrac');
            end
            if isscalar(t.sFrac)
                sFrac=t.sFrac; %#ok<PROP>
                t.sFrac = zeros(size(t.neurIdentities));
                t.sFrac(:) = sFrac; %#ok<PROP>
            end
            
            this=t;
            
        end
        % ---------------------------------------------------------------
        function [tauM,tauD,tauS,p0,rMax,Iwidth,Ith] = returnOutputs(this)
            if this.useGPU
                tauM=gpuArray(single(this.tauM)); 
                tauD=gpuArray(single(this.tauD)); 
                tauS=gpuArray(single(this.tauS));
                p0=gpuArray(single(this.p0)); 
                rMax=gpuArray(single(this.rMax));
                Iwidth=gpuArray(single(this.Iwidth));
                Ith=gpuArray(single(this.Ith));
            else
                tauM=this.tauM; 
                tauD=this.tauD; 
                tauS=this.tauS;
                p0=this.p0; 
                rMax=this.rMax;
                Iwidth=this.Iwidth;
                Ith=this.Ith;
            end
        end
        % ---------------------------------------------------------------
    end
    
    % OUTPUT PROPERTIES
    properties (SetAccess = protected)
        % Time Constants
       tauM;
       tauD;
       tauS;
       
       % synaptic fraction
       p0;
       
       % Max Firing rates
       rMax;
       
       % Current response
       Ith;
       Iwidth;
    end
    
end