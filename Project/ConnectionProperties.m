classdef ConnectionProperties < ConnectionInterface
    % Implementation of connections class, that determines weights for our
    % neural units. A different implementation is called for 2D topography
    % of cells. This connection code can support assigning connection to an
    % arbitrayr number of connection types, not just EE, EI, and so on. If
    % you place new types into the types property, you can then control
    % connections between that many populations.  I wrote this to be
    % general enough such that if I had an idea for changing the way I
    % connect things, I wouldn't have to reorganize my code much. And I
    % prefer the order that object-oriented code brings to the table.
    
    %% INPUT PROPERTIMES
    % Specific to this implementation of the the connections interface
    % class.
    properties (Access = public)  
        % Parameters that specify the asymmetric and reccurent connection
        % properties, and as well as standard ...
        
        % ----------------------------------------------------------------
        % OPTION 1
            % Base connectivity level for the particular type
            Base;
            % Base.EE = base EE value before randomness componenet
            % Base.IE = base IE value before randomness componenet
            % et cetera
       % ----------------------------------------------------------------% ----------------------------------------------------------------
        % OPTION 2
            % A connection type either receives an asymmetric input and
            % recurrent input for it's population .. OR ... it gets a base
            Asym;
            Rec;
            % Asym.EE = asymmetic EE value before randomness componenet
            % Asym.IE = asymmetric IE value before randomness componenet
            % et cetera
        % ----------------------------------------------------------------
        
        % A struct whose members specify the randomness in the connection
        % types above!
        Sigma;
        
        % Recurrence Mode - controls whether the recurrence is applied onto
        % the unit's self, or whether there are a population of units
        % defined as a group, who then are recurrent onto themselves
        PopulationMode;
        fracPop;
        
        % Can veto GPU mode
        % (Output vectors are stored on the GPU so that downstream
        % simulation can take advantage of it.)
        disableGPU=false;
    end
    properties (Access = private)
        % Allowable population types
        types = {'EE','EI','II','IE'};
        
        % Vectors that track exc and inhibitory locations
        exc; inh;
    end
    
    %% METHODS 
    methods
        % ------------------------------------------------------------
        function [this]=ConnectionProperties(this,varargin) %#ok<INUSL>
            % Trackers
            p_set = false;
            % Apply optional inputs
            for v = 1:2:numel(varargin)
                switch varargin{v}
                    case 'Base', this.Base=varargin{v+1};
                    case 'Asym', this.Asym=varargin{v+1};
                    case 'Rec', this.Rec=varargin{v+1};
                    case 'Sigma', this.Sigma=varargin{v+1};
                    case 'PopulationMode' 
                        p_set = true;
                        this.PopulationMode=varargin{v+1};
                    otherwise, error('Unrecognized input');
                end
            end
            
            % If certain inputs were set or not, we apply some basic checks
            % and balances
            if ~p_set
                for t = this.types
                    this.PopulationMode.(t{1}) = false;
                end
            end
            this.checkTypes();
        end
        % ------------------------------------------------------------
        function [W] = returnOutputs(this)
            % Common to all interfaces, there is a return function that
            % returns all vectors that will be used in the simulation. Some
            % objects return one output; some many.
            W=this.W;
        end
        % ------------------------------------------------------------
        function [this]=generateConnections(this)
            % This is the function that calculates the output parameters
            % that the downstream simulation uses.

            % Check that inputted settings look fine
            this.checkTypes();
            % Setup shortcut alias
            t=this;
            
            % Create random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');
            
            % Get indices of the exitatory and inhibitory cells
            this.exc         = t.neurIdentities == 0;
            this.inh         = t.neurIdentities == 1;
            nNeurons    = numel(t.neurIdentities);
            
            % Initialize all of the various weight matrices to be the size
            % of the eventual W matrix
            W.EE = zeros(nNeurons,nNeurons);
            W.EI = zeros(nNeurons,nNeurons);
            W.II = zeros(nNeurons,nNeurons);
            W.IE = zeros(nNeurons,nNeurons);
            % Then the eventual recepticle of all these instatiated values
            t.W = zeros(nNeurons,nNeurons);
            
            %% First, compute any that have base probabilities set
            for f = fields(this.Base)
                type=f{1};
                sigma=0;
                if isfield(this.Sigma,type)
                    sigma = this.Sigma.(type);
                    if ischar(sigma)
                        % If it's a string, the user inputted a function to
                        % determine the sigma value fo the function.
                        % sigma=this.conditionalParam(sigma);
                    elseif isstruct(sigma)
                        % If it evaluates to a struct, that's a sign that
                        % the user has provided information about Rec or
                        % Asym for this particular connection type, and so
                        % we should skip it
                        continue;
                    end
                end
                
                [popX,popY]=this.whichElements(type);
                base = this.Base.(type);
                if isnumeric(base)
                    W.(type)(popX,popY) = base + ...
                        sigma*rand(r,sum(popX),sum(popY));
                end
            end
            
            %% Higher precedence settings, are the asymmetric/recurrent modes
            % which do not apply blanket probabilities to the entire population
            % type.
            if ~isempty(this.Asym)
                for f = fields(this.Asym)
                    type=f{1};
                    [popX,popY]=this.whichElements(type);

                    % Get the asymmetric sigma
                    if isfield(this.Sigma,'Asym') && ...
                            isfield(this.Sigma.Asym, type)
                        sigmaAsym = this.Sigma.Asym.(type);
                    else
                        sigmaAsym = 0;
                    end
                    % Get the recurrent sigma
                    if isfield(this.Sigma,'Rec') && ...
                            isfield(this.Sigma.Rec, type)
                        sigmaRec = this.Sigma.Rec.(type);
                    else
                        sigmaRec = 0;
                    end
                    % Get the asymmetric base
                    baseAsym = this.Asym.(type);
                    % Get the recurrent base
                    baseRec = this.Rec.(type);

                    % Check if pop mode set
                    popmode = this.PopulationMode(type);
                    if popmode
                        % To fill
                    else % Single unit recurrence mode
                        W.(type)(popX,popY) = baseAsym + ...
                            sigmaAsym*rand(r,sum(popX),sum(popY));
                        recInd = ...
                            ConnectionProperties.getDiagonal(popX);
                        W.(type)(recInd) = baseRec + ...
                            sigmaRec*sigmaAsym*rand(r,1,numel(recInd));
                    end  
                end
            end
            
            %% Set total connection matrix
            try
                if gpuDeviceCount > 0
                    t.W = gpuArray(W.EE + W.IE + W.II + W.EI);
                else
                    t.W = gpuArray(W.EE + W.IE + W.II + W.EI);
                end
            catch
                % if gpuDeviceCount function not found or if gpuArray fails
                % to initialize because of CUDA driver problems, then it
                % catches here and simply creates a standard array.
                warning('Caught error while trying to assign output');
                t.W = gpuArray(W.EE + W.IE + W.II + W.EI);
            end
            
            % Set the matrix that can be optionally returned
            W=t.W; % for return purposes only
            this=t;
        end
        % ------------------------------------------------------------
        function checkTypes(this)
            % Checks that inputs that use has inputted sensible types
            checks = {'Base','Rec','Asym'};
            for c = checks 
                if ~isempty(this.(c{1}))
                    if sum(~ismember(fields(this.(c{1})),this.types))
                        error('Wrong field inputs to %s',c{1});
                    end
                end
                % If the user provided any functions for the parameters
                % instead of hard coding them to a property -- then they're
                % applied here
                if ~isempty(this.(c{1}))
                    for f = fields(this.(c{1}))
                        if ischar(this.(c{1}).(f{1}))
                            this.(c{1}).(f{1}) = ...
                                this.conditionalParam(this.(c{1}).(f{1}));
                        end
                    end
                end
            end
            
            % Ensure that every field given to Asym is also provided to Rec
            if ~isempty(this.Asym)
                for f = fields(this.Asym)
                    if ~isfield(this.Rec,f{1})
                        error('Field %s set in Asym not set in Rec',f{1});
                    end
                end
            end
            
            % Vice versa
            if ~isempty(this.Asym)
                for f = fields(this.Rec)
                    if ~isfield(this.Asym,f{1})
                        error('Field %s set in Rec not set in Asym',f{1});
                    end
                end
            end
            
        end
        % ------------------------------------------------------------
        function out = conditionalParam(this,str) %#ok<INUSL>
            out = 0;
            try
                eval(['out = ' str]);
            catch
                error('The functional input was incorrect.');
            end
        end
        % ------------------------------------------------------------
        function [popX, popY] = whichElements(this,type)
            loc = find(ismember(type,this.types));
            switch loc
                case 1
                    popX=this.exc; popY=this.exc;
                case 2
                    popX=this.exc; popY=this.inh;
                case 3
                    popX=this.inh; popY=this.inh;
                case 4
                    popX=this.inh; popY=this.exc;
                otherwise
                    error('Invalid option');
            end
        end
    end
    % ------------------------------------------------------------
    methods (Static)
        function out = getDiagonal(elements)
             % Find the diagonal indices for a subset of cells
            out = sub2ind(size(t.W),find(elements),find(elements));
        end
    end
    
end