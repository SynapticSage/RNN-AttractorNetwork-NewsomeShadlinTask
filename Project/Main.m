% Main for the script that implements an attractor based neural firing
% model on a stimulus regime that's supposed to mimic the Mante task.

%% General Paramters
% The lines below this comment provide an optional entry point for my
% ParameterExplorer class to take over execution of the script, else, it
% simply runs the script with the following default parameters. For reasons
% why using this parameter explorer is smart, for one, it automatically
% stores data in a folder system structured like the parameter space
% explored. Second, it is able to run many instances of this script in
% parallel, as many as the processor can support. Third, it has an ability
% to stop/start in the middle of paramter exploration (checkpoints,
% automatically).
if ~exist('params','var') || ~isstruct(params) 
   params = struct( ...
           ... Model Resolution/Precision
           'dt',    2e-4, ... Model bin size
           'bin_t', 5e-2, ... How to bin post hoc for analysis
       ... ---TASK PARAMETERS---
           ... General Trial Controls
           'nTrials',       100, ...
           'trialDuration', 3, ...
           ... Context Stimulus
           'context_trialStart',    0.3, ...
           'context_trialEnd',      1.1, ...
           ... Color Motion Stimulus
           'color_trialStart',      0.35, ...
           'color_trialEnd',      1.1, ...
           ... Dot Color Stimulus
           'dot_trialStart',        0.35, ...
           'dot_trialEnd',      1.1, ...
       ... ---NEURAL NETWORK PARAMETERS---
           ... Number of neurons of each type
           'nInh',      5, ...
           'nExc',      40, ...
           ... Neural Properties
           'rMax0E',    100, ...
           'rMax0I',    200, ...
           'p0E',       0.3, ...
           'p0I',       0.1, ...
           ...Time Constants
           'tausE',     0.025, ...
           'tausI',     0.005, ...
           'tauDbar',   0.025, ...
           'tauDvar',   0,     ...
           'taum',      0.010, ...
           ... Connection Properties
           'Wrecurrent',    200, ...
           'sigmaWEE',      0, ...
           'sigmaWasym',    0, ...
           'WIEval',        -320, ...
           'sigmaIE',       0 ...
       );
end

% Temporary note - properties that have to be defined post-parameter stage
% not yet accounted for: (Ith_i, Ith_e, IsigmaTh), (Iwidth_i Iwidth_e),
% WEEasym WEIval pmean

% Setup neuron counts
neurIdentites=[false(1,params.nExc) true(1,params.nInh)];

%% Initialize Parameter-based Simulation Vectors
% ------------------------------------------------------------------------
% NEURAL VECTORS
N = NeuronProperties; % Encapsulated code for setting up overall neuron/synaptic vectors
    % First, we set the properties to build the the output properties
    N.rMax0E=params.rMax0E; N.rMax0I=params.rMax0I;
    N.p0E=params.p0E; N.p0I=params.p0I;
    N.tausE=params.tausE; N.tausI=params.tausI;
    N.tauDbar=params.tauDbar; N.tauDvar=params.tauDvar;
    % Then, we invoke the generation method to create them
    N.generateOutputParams;

% ------------------------------------------------------------------------
% CONECTION VECTORS
C = ConnectionProperties; % Encapsulated code for computing overall W vector
    % First, we set the properties to build the W matrix
    C.Rec.EE = params.Wrecurrent;
    C.Sig.Rec.EE = params.sigmaWEE; C.Sig.Asym.EE = params.sigmaWEE;
    C.Base.IE = params.WIEval; C.Sigma.IE = params.sigmaIE;
    % then, we invoke the generation method to create W
    W = C.generateConnections();

% ------------------------------------------------------------------------
% STIMULI VECTORS
    % (1) Setup Context Stimulus
    % (2) Setup Dot Stimulus
    % (3) Setup Color Stimulus

%% Execute Simulation
for trial = 1:max_trials
    trial
    r = zeros(length(t),Ncells);    % Firing rate for each cell at all time points
    D = zeros(length(t),Ncells);    % Depression variable for each cell at all time points
    S = zeros(length(t),Ncells);    % Synaptic gating variable for each cell at all time points
    
    %% Per Stimulation
    for stim = 1:Nstims
        
        if ( trial_reset || stim == 0 ) || ...
                ( multistim && mod(stim,Nmax+1)==0 )
            r(1+Nt*(stim-1),:) = 0.0;                            % Initializing if resetting to different stimuli
            D(1+Nt*(stim-1),:) = 1.0;
            S(1+Nt*(stim-1),:) = 0.0;
        else
            r(1+Nt*(stim-1),:) = r(Nt*(stim),:) ;               % Do not initialize if continuing to count stimuli
            D(1+Nt*(stim-1),:) = D(Nt*(stim),:);
            S(1+Nt*(stim-1),:) = S(Nt*(stim),:);
        end
        
        %% Step Through Times
        for i = 2+Nt*(stim-1):Nt*(stim)                            % Now integrate through time
            I = S(i-1,:)*W+Iapp(i,:) ...        % I depends on feedback (W*S) and applied current
                + sigma*randn(s1,1)/sqrt(dt);      % and additional noise
            % S(:,i-1) is the vector of synaptic gating
            % from the previous time step for all cells
            % This gets multiplied by the weight matrix
            % to give total feedback current.
            
            rinf = rmax./(1.+exp(-(I-Ith)./Iwidth));        % Firing rate curve gives the steady state r
            r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/taum);  % Update r from the previous timestep
            
            Dinf = 1./(1.+p0.*r(i-1,:).*taud);                  % Steady state value of D for Poisson spiking
            D(i,:) = Dinf + ( D(i-1,:)-Dinf).*...
                exp(-dt*(p0.*r(i-1,:)+1./taud));  % Update with adjusted time constant
            
            Sinf = sfrac*p0.*r(i,:).*D(i,:).*taus./(1.0+sfrac*p0.*r(i,:).*D(i,:).*taus); % Steady state value of synaptic gating vatiable assuming vesicle release at a rate p0*r*D
            S(i,:) = Sinf + ( S(i-1,:)-Sinf).*...
                exp(-dt*(sfrac*p0.*r(i,:).*D(i,:)+1./taus)); % update S with adjusted tau
        end % continue to next time step

        if (multistim == false  && stim > 0 ) || mod(stim,Nmax+1) > 0 
            
            % First half of trials obtain mean network responses to
            % stimuli, used later for confusibility matrix
            if trial <= Num_of_trials
                meanresponse1(stim,:) = meanresponse1(stim,:) + ...
                    mean(r(i-Nsec1:i,:));
                sdresponse1(stim,:) = sdresponse1(stim,:) + ...
                    mean(r(i-Nsec1:i,:)).*mean(r(i-Nsec1:i,:));
            else
                % Second half of trials used as test responses
                mresponse1(stim,:,trial-Num_of_trials) =  mean(r(i-Nsec1:i,:));
            end
        end
        
    end 
    
    %% Post-trial plotting
    if figureson
        figure(1)
        imagesc(r(:,1:end-1)')
        colorbar
        drawnow
    end
    
    meanrate = meanrate + r;
    stdrate = stdrate + r.*r;
    
end

%% Post-simulation Analysis