classdef NeuronProperties < UnitsInterface
    
    properties
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
        tauDvar;
        
        % Specification for fraction of receptors activated by release
        sfrac=1;
    end
    
    properties (Access=protected)
        initFunction; % Stores the initialization function used to initialize output variables. This controls whether or not it creates gpu arrays
    end
    
    methods
        % ---------------------------------------------------------------
        function this = NeuronProperties(this)
            % Constructor method that primarily determine how the generator
            % should initialize it's data types, gpu-driven or not
            
            fprintf('Set the neuron properties ...\n');
            
            try
            if gpuDeviceCount > 0
               this.initFunction = @(s) gpuArray(zeros(size(s)));
            else
               this.initFunction = @(s) zeros(size(s));
            end
            catch
                % if gpuDeviceCount function not found or if gpuArray fails
                % to initialize because of cuda driver problems, then it
                % catches here and simply creates a standard array.
                this.initFunction = @(s) zeros(size(s));
            end
            
        end
        % ---------------------------------------------------------------
        function this= generateConnections(this)
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
            exc         = t.identities == 0;
            inh         = t.identities == 1;
            
            % Assign the rMax vector
            t.rMax = this.initFunction(t.neurIdentities);
            t.rMax(exc) = t.rMax0E;
            t.rMax(inh) = t.rMax0I;
            
            % Assign the synaptic time const
            t.tauS = this.initFunction(t.neurIdentities);
            t.tauS(exc) = t.tausE;
            t.tauS(inh) = t.tausI;
            
            % Assign the depression time const...
            % creates a random constant for the excitatory cells, and a
            % non-random for inhibitory cells
            t.tauD = zeros(size(t.neurIdentities));
            t.tauD(exc) = t.tauDbar + t.tauDvar*rand(r,size(t.tauD(exc)));
            t.tauD(inh) = t.tauDbar;
            
            % Specify the release propabilities per synapse
            t.p0 = this.initFunction(t.neurIdentities);
            t.p0(exc) = t.p0E;
            t.p0(inh) = t.p0I;
            
            % Specify the sfrac
            if isscalar(t.sFrac)
                sFrac=t.sFrac;
                t.sFrac = this.initFunction(t.neurIdentities);
                t.sFrac(:) = sFrac;
            end
            
        end
    end
    
end