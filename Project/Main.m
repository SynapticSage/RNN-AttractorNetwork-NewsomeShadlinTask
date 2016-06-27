% Main for the script that implements an attractor based neural firing
% model on a stimulus regime that's supposed to mimic the Mante task.

%% General Flags
useGPU = false;
trialReset = false;

%% Pre-processing
if useGPU
    reset(gpuDevice(1)); % Resets gpu memory, if anything is in it.
end

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
           'dt',    2e-3, ... Model bin size
           'bin_t', 5e-2, ... How to bin post hoc for analysis
       ... ---TASK PARAMETERS---
           ... General Trial Controls
           'nTrials',       5, ...
           'trialDuration', 3, ...
           ... General Stimulus
           'iFrac',         0.33, ...   Randomly select a third of the population for a stimulus to receive
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
           'tauM',      0.010, ...
           ... Input Properties
           'IwidthE',   3, ...
           'IthE',      18, ...
           'IwidthI',   5, ...
           'IthI',      20, ...
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
nCells = params.nExc + params.nInh;

%% Initialize Parameter-based Simulation Vectors
% ------------------------------------------------------------------------
% NEURAL VECTORS
NP = NeuronProperties; % Encapsulated code for setting up overall neuron/synaptic vectors
    % First, we set the properties to build the the output properties
    NP.neurIdentities = neurIdentites;
    NP.rMax0E=params.rMax0E; NP.rMax0I=params.rMax0I;
    NP.p0E=params.p0E; NP.p0I=params.p0I;
    NP.tausE=params.tausE; NP.tausI=params.tausI;
    NP.tauDbar=params.tauDbar; NP.tauDvar=params.tauDvar;
    NP.tauMbar=params.tauM;
    NP.IwidthI =  params.IwidthI; NP.IwidthE = params.IwidthE;
    NP.Ith0E =  params.IthE; NP.Ith0I = params.IthI;

    % Then, we invoke the generation method to create them
    NP = NP.generateOutputParams;

    % last return the gpu vectors for the simulation
    [tauM,tauD,tauS,p0,rMax,Iwidth,Ith] = NP.returnOutputs(); clear NP;

% ------------------------------------------------------------------------
% CONECTION VECTORS
Connect = ConnectionProperties; % Encapsulated code for computing overall W vector
    % First, we set the properties to build the W matrix
    Connect.neurIdentities=neurIdentites;
    Connect.Rec.EE = params.Wrecurrent;
    Connect.Sigma.Rec.EE = params.sigmaWEE;
    Connect.Sigma.Asym.EE = params.sigmaWEE;
    Connect.Base.IE = params.WIEval; Connect.Sigma.IE = params.sigmaIE;

    % then, we invoke the generation method to create W
    Connect = Connect.generateConnections();

    % last return the gpu vectors for the simulation
    [W] = Connect.returnOutputs(); clear Connect;

% ------------------------------------------------------------------------
% STIMULI VECTORS
    GeneralStim = InputStimulus_Simple();
        % The following parameters should be put into the moving parameter
        % set
        GeneralStim.neurIdentities = neurIdentites;
        GeneralStim.dt = params.dt;
        GeneralStim.nTrials = params.nTrials; % number of trials
        GeneralStim.trialDuration = params.trialDuration; % seconds
        GeneralStim.iFrac = 0.3;
    % (1) Setup Context Stimulus
    Context = GeneralStim;
        Context.trialStart  = params.context_trialStart;
        Context.trialStop   = params.context_trialEnd;
        Context.nStates     = 2;
        Context=Context.generateStimuli();
        [Iapp_context, trialBin, t] = Context.returnOutputs();
    % (2) Setup Dot Stimulus
    Dots = GeneralStim;
        Dots.trialStart     = params.dot_trialStart;
        Dots.trialStop      = params.dot_trialEnd;
        Dots.nStates        = 3;
        Dots=Dots.generateStimuli();
        [Iapp_dots] = Dots.returnOutputs();
    % (3) Setup Color Stimulus
    Color = GeneralStim;
        Color.trialStart    = params.color_trialStart;
        Color.trialStop     = params.color_trialEnd;
        Color.nStates        = 3;
        Color=Color.generateStimuli();
        [Iapp_color] = Color.returnOutputs();
    clear Context Dots Color GeneralStim;
    clear Context Dots Color GeneralStim;

   % Store the superposition of the inputs;
   Iapp = Iapp_context + Iapp_color + Iapp_dots;

%% Initialize statistic trackers
% meanrate        = zeros(length(t),nCells);                     % Mean time-dependence averaged across trials
% stdrate         = zeros(length(t),nCells);                      % Std of time-dependence across trials
% mresponse1      = zeros(Nstims,nCells,Num_of_trials);        % Response after time-averaging to each stimulus
% meanresponse1   = zeros(Nstims,nCells);                   % Mean response after averaging across trials
% sdresponse1     = zeros(Nstims,nCells);                     % Std of responses after averaging across trials

%% Execute Simulation
for trial = 1:params.nTrials

    r = zeros(length(t),nCells);    % Firing rate for each cell at all time points
    D = zeros(length(t),nCells);    % Depression variable for each cell at all time points
    S = zeros(length(t),nCells);    % Synaptic gating variable for each cell at all time points

    if useGPU
        r=gpuArray(r);
        D=gpuArray(D);
        S=gpuArray(S);
    end

         if trialReset
            r(1+Nt*(stim-1),:) = 0.0;                            % Initializing if resetting to different stimuli
            D(1+Nt*(stim-1),:) = 1.0;
            S(1+Nt*(stim-1),:) = 0.0;
%         else
%             r(1+Nt*(stim-1),:) = r(Nt*(stim),:) ;               % Do not initialize if continuing to count stimuli
%             D(1+Nt*(stim-1),:) = D(Nt*(stim),:);
%             S(1+Nt*(stim-1),:) = S(Nt*(stim),:);
        end

    %% Step Through Times
    startInd   = trialBin(trial,1);
    stopInd    = trialBin(trial,2);
    for i = startInd:stopInd

        I = S(i-1,:)*W+Iapp(i,:) ...        % I depends on feedback (W*S) and applied current
            + sigma*randn(s1,1)/sqrt(dt);   % and additional noise

        rinf = rmax./(1.+exp(-(I-Ith)./Iwidth));        % Firing rate curve gives the steady state r
        r(i,:) = rinf + (r(i-1,:)-rinf)*exp(-dt/tauM);  % Update r from the previous timestep

        Dinf = 1./(1.+p0.*r(i-1,:).*tauD);                  % Steady state value of D for Poisson spiking
        D(i,:) = Dinf + ( D(i-1,:)-Dinf).*...
            exp(-dt*(p0.*r(i-1,:)+1./tauD));  % Update with adjusted time constant

        Sinf = sfrac*p0.*r(i,:).*D(i,:).*tauS./(1.0+sfrac*p0.*r(i,:).*D(i,:).*taus); % Steady state value of synaptic gating vatiable assuming vesicle release at a rate p0*r*D
        S(i,:) = Sinf + ( S(i-1,:)-Sinf).*...
            exp(-dt*(sfrac*p0.*r(i,:).*D(i,:)+1./taus)); % update S with adjusted tau
    end

    %% Add to post-trial statistics
%     if (multistim == false  && stim > 0 ) || mod(stim,Nmax+1) > 0
%
%         % First half of trials obtain mean network responses to
%         % stimuli, used later for confusibility matrix
%         if trial <= Num_of_trials
%             meanresponse1(stim,:) = meanresponse1(stim,:) + ...
%                 mean(r(i-Nsec1:i,:));
%             sdresponse1(stim,:) = sdresponse1(stim,:) + ...
%                 mean(r(i-Nsec1:i,:)).*mean(r(i-Nsec1:i,:));
%         else
%             % Second half of trials used as test responses
%             mresponse1(stim,:,trial-Num_of_trials) =  mean(r(i-Nsec1:i,:));
%         end
%     end

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
