classdef ConnectionProperties
    % Implementation of connections class, that determines weights for our
    % neural units. A different implementation is called for 2D topography
    % of cells. This connection code can support assigning connection to an
    % arbitrayr number of connection types, not just EE, EI, and so on. If
    % you place new types into the types property, you can then control
    % connections between that many populations.  I wrote this to be
    % general enough such that if I had an idea for changing the way I
    % connect things, I wouldn't have to reorganize my code much. And I
    % prefer the order that object-oriented code brings to the table.

    %% INPUT Variables
    % Specific to this implementation of the the connections interface
    % class.
    properties (Access = public)
        neurIdentities; % Stores information about the neural units/groups, how they present across the vector -- it's how we label which neurons are which types, and that determines by the rules in the assignment method how to setup the rules.

        % Parameters that specify the asymmetric and reccurent connection
        % properties, and as well as standard ...
        % ----------------------------------------------------------------
        % OPTION 1
            % Base connectivity level for the particular type
            Base    = struct('EE',0,'IE',0,'EI',0,'II',0)
            % Struct tracking probability of connections between populations
            P       = struct('EE',1,'IE',1,'EI',1,'II',1);
       % ----------------------------------------------------------------% ----------------------------------------------------------------
        % OPTION 2
            % A connection type either receives an asymmetric input and
            % recurrent input for it's population .. OR ... it gets a base
            Rec     = struct('EE',0,'II',0)
        % ----------------------------------------------------------------

        % A struct whose members specify the randomness in the connection
        % types above!
        Sigma = struct('EE',0,'IE',0,'EI',0,'II',0);

        % Recurrence Mode - controls whether the recurrence is applied onto
        % the unit's self, or whether there are a population of units
        % defined as a group, who then are recurrent onto themselves .. if
        % these fields are set to a value other than false, then they
        % function as probabilities that control what's defined as a
        % recurrent population
        PopulationMode=struct('EE',false,'II',false);

        % Can veto GPU mode
        % (Output vectors are stored on the GPU so that downstream
        % simulation can take advantage of it.)
        useGPU=false;
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
        function [this]=ConnectionProperties(~,varargin)

            % If certain inputs were set or not, we apply some basic checks
            % and balances
            for t = this.types
                this.PopulationMode.(t{1}) = false;
            end

            this.checkTypes();
        end
        % ------------------------------------------------------------
        function [W] = returnOutputs(this)
            % Common to all interfaces, there is a return function that
            % returns all vectors that will be used in the simulation. Some
            % objects return one output; some many.

            if this.useGPU
                W=gpuArray(single(this.W));
            else
                W=this.W;
            end
        end
        % ------------------------------------------------------------
        function [this]=generateConnections(this)
            % This is the function that calculates the output parameters
            % that the downstream simulation uses.

            % Check that inputted settings look fine
            this.checkTypes();

            % Create random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');

            % Get indices of the exitatory and inhibitory cells
            this.exc         = this.neurIdentities == 0;
            this.inh         = this.neurIdentities == 1;
            nNeurons         = numel(this.neurIdentities);

            % Initialize all of the various weight matrices to be the size
            % of the eventual W matrix
            W.EE = zeros(nNeurons,nNeurons); %#ok<*PROP>
            W.EI = zeros(nNeurons,nNeurons);
            W.II = zeros(nNeurons,nNeurons);
            W.IE = zeros(nNeurons,nNeurons);
            % Then the eventual recepticle of all these instatiated values
            this.W = zeros(nNeurons,nNeurons);

            %% First, compute any that have base probabilities set
            for f = fields(this.Base)'
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

                [popIndices,cScale]=this.return_XtoY(type);

                base = this.Base.(type);
                if isnumeric(base)
                    W.(type)(popIndices) = cScale*(...
                        base + ...
                        sigma*rand(r,size(W.(type)(popIndices)))...
                        );
                end
            end

            %% Higher precedence settings, are the asymmetric/recurrent modes
            % which do not apply blanket probabilities to the entire population
            % type.
            if ~isempty(this.Rec)
                for f = fields(this.Rec)'
                    % connection type string, eg 'EI'
                    type=f{1};
                    % Sigma Reccurrence
                    if isfield(this.Sigma,'Rec') && ...
                            isfield(this.Sigma.Rec, type)
                        sigmaRec = this.Sigma.Rec.(type);
                    else
                        sigmaRec = 0;
                    end
                    % Recurrence Base Value
                    baseRec = this.Rec.(type);
                    % Populate Mode?
                    % Check if pop mode set
                    popmode = this.PopulationMode.(type);
                    if popmode
                        % To fill
                    else % Single unit recurrence mode
                        [~,~,popX,~]=this.return_XtoY(type);
                        recInd = ...
                            ConnectionProperties.getDiagonal(W.(type),popX);
                        W.(type)(recInd) = baseRec + ...
                            sigmaRec*rand(r,1,numel(recInd));
                    end
                end
            end

            %% Set total connection matrix
            this.W = W.EE + W.IE + W.II + W.EI;

        end
        % ------------------------------------------------------------
        function [popIndices, averageOut, popX, popY] = ...
                return_XtoY(this,type)

            % Obtain logical vectors specifying the from (X) and to (Y)
            % populations -- ie from neurons labeld by popX to neurons
            % labeled by popY
            loc = find(ismember(this.types,type));
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

            % generate matrix of 1's that detail who is about to connect to
            % who
            popMat = zeros(numel(popX));

            popMat(popX,popY) = true;
            popOnCopy = popMat;

            if this.P.(type)
              % randomly prune connection according to any probabilities
              % given
              randMat   = rand(numel(popX));
              popMat  = popMat .* randMat;
              popMat  = popMat <= this.P.(type);

              % Finalize component for these randomized connections
              popMat = logical(popMat .* popOnCopy);

              % Last compute the cScale factor by the average number of
              % connections extended out by a given cell
              averageOut = sum(popMat,2);
              averageOut = 1/mean(averageOut);
            else

            end

            % now we need to calculate the indices that follow as a result
            % of the population X projecting to population Y.
            popIndices = reshape(popMat,1,[]);
        end
        % ------------------------------------------------------------
        function out = visualize(this,savedir)

            fprintf('Value of this.W:\n');
            disp(this.W);

            % Prints a visual description of the generated connections into
            % a save folder
            if nargin==1;savedir=pwd;end;

            % Generate plot of connections
            out=figure(100); clf;
            subplot(1,2,1);
            axs = 1:numel(this.neurIdentities);
            laxs=[axs(1)-0.5 axs axs(end)+0.5];
            imagesc(axs,axs,this.W);
            colorbar;
            title('Connectivity Matrix');
            xlabel('Neuron #');
            ylabel('Neuron #');
            % Demarcate EI border, if neurIdentites were grouped by E and
            % I.
            loc = diff(this.neurIdentities ~=0);
            if sum(loc) == 1
                loc = find(loc)+0.5;
                hold on;
                l=line(repmat(loc,size(laxs)),laxs);
                set(l,'color','r','linewidth',2,'linestyle',':');
                l=line(laxs,repmat(loc,size(laxs)));
                set(l,'color','r','linewidth',2,'linestyle',':');
            end

            subplot(1,2,2);
            laxs = 1:numel(this.neurIdentities);
            logmeasure=log10(abs(this.W));
            logmeasure(logmeasure==-inf) = 0;
            logmeasure = logmeasure.*sign(this.W);
            imagesc(axs,axs,logmeasure);
            colorbar;
            title('''SignedLog'' Connectivity Measure');
            xlabel('Neuron #');
            ylabel('Neuron #');
            % Demarcate EI border, if neurIdentites were grouped by E and
            % I.
            loc = diff(this.neurIdentities) ~=0;
            if sum(loc) == 1
                loc = find(loc)+0.5;
                hold on;
                l=line(repmat(loc,size(laxs)),laxs);
                set(l,'color','r','linewidth',2,'linestyle',':');
                l=line(laxs,repmat(loc,size(laxs)));
                set(l,'color','r','linewidth',2,'linestyle',':');
            end

            % Save
            warning off;
            saveThis(out,savedir,'Connections','fig');
            saveThis(out,savedir,'Connections','png');
            warning on;

        end
        % ------------------------------------------------------------
        function checkTypes(this)
            % Checks that inputs that use has inputted sensible types
            checks = {'Base','Rec'};
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
                    for f = fields(this.(c{1}))'
                        if ischar(this.(c{1}).(f{1}))
                            this.(c{1}).(f{1}) = ...
                                this.conditionalParam(this.(c{1}).(f{1}));
                        end
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
    end
    % ------------------------------------------------------------
    methods (Static)
        function out = getDiagonal(W,elements)
             % Find the diagonal indices for a subset of cells
            out = sub2ind(size(W),find(elements),find(elements));
        end
    end

    %% OUTPUT Variables
    properties (SetAccess = private)
        W; % Stores the resulting connection strenghts
    end

end
