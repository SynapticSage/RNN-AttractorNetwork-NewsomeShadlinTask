% Main for the script that implements an attractor based neural firing
% model on a stimulus regime that's supposed to mimic the Mante task.

%% General Flags
useGPU       = false;
trialReset   = false;
figureson    = true;
PE_on        = true;

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
if exist('params','var') && parameterexplorer_on
  % Initialize parameters encoded in params struct
  ParameterExplorer.swapParamSet(params);
  % Automatically set save location based on parameters
  savedir = ParameterExplorer.savelocation(params,...
      'projectfolder',projectfolder);
  fprintf('SaveLocation: %s\n',savedir);
else
  savedir = '~/Data/Miller/Miscelleaneous    '
  Load_DefaultParams;
end

% Setup neuron counts
neurIdentites=[false(1,nExc) true(1,nInh)];
nCells = nExc + nInh;

% Setup random stream for main script (different than random streams that
% are operating to construct stimuli, connections, and neuron properties,
% respectively)
MainStream = RandStream.create('mrg32k3a','NumStreams',1,'seed','shuffle');

%% Initialize Parameter-based Simulation Vectors
% ------------------------------------------------------------------------
% NEURAL VECTORS
NP = NeuronProperties; % Encapsulated code for setting up overall neuron/synaptic vectors
    % First, we set the properties to build the the output properties
    NP.neurIdentities = neurIdentites;
    NP.rMax0E=rMax0E; NP.rMax0I= rMax0I;
    NP.p0E= p0E; NP.p0I=p0I;
    NP.tausE=tausE; NP.tausI=tausI;
    NP.tauDbar=tauDbar; NP.tauDvar=tauDvar;
    NP.tauMbar=tauM;
    NP.IwidthI =  IwidthI; NP.IwidthE = IwidthE;
    NP.Ith0E =  IthE; NP.Ith0I = IthI;

    % Then, we invoke the generation method to create them
    NP = NP.generateOutputParams;

    % last return the gpu vectors for the simulation
    [tauM,tauD,tauS,p0,sFrac,rMax,Iwidth,Ith] = NP.returnOutputs();
    clear NP;

% ------------------------------------------------------------------------
% CONECTION VECTORS
Connect = ConnectionProperties; % Encapsulated code for computing overall W vector
    % First, we set the properties to build the W matrix
    Connect.neurIdentities=neurIdentites;
    Connect.Rec.EE = Wrecurrent;
    Connect.Sigma.Rec.EE = sigmaWEE;
    Connect.Sigma.Asym.EE = sigmaWEE;
    Connect.Base.IE = WIEval; Connect.Sigma.IE = sigmaIE;

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
        GeneralStim.dt = dt;
        GeneralStim.nTrials = nTrials; % number of trials
        GeneralStim.trialDuration = trialDuration; % seconds
        GeneralStim.iFrac = 0.3;
    % (1) Setup Context Stimulus
    Context = GeneralStim;
        Context.trialStart  = context_trialStart;
        Context.trialStop   = context_trialEnd;
        Context.nStates     = 2;
        Context=Context.generateStimuli();
        [Iapp_context, trialBin, t] = Context.returnOutputs();
    % (2) Setup Dot Stimulus
    Dots = GeneralStim;
        Dots.trialStart     = dot_trialStart;
        Dots.trialStop      = dot_trialEnd;
        Dots.nStates        = 3;
        Dots=Dots.generateStimuli();
        [Iapp_dots] = Dots.returnOutputs();
    % (3) Setup Color Stimulus
    Color = GeneralStim;
        Color.trialStart    = color_trialStart;
        Color.trialStop     = color_trialEnd;
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

r = zeros(length(t),nCells);    % Firing rate for each cell at all time points
D = zeros(length(t),nCells);    % Depression variable for each cell at all time points
S = zeros(length(t),nCells);    % Synaptic gating variable for each cell at all time points
if useGPU
        r=gpuArray(r);
        D=gpuArray(D);
        S=gpuArray(S);
end
for trial = 1:nTrials

     if trialReset
        r(1+Nt*(stim-1),:) = 0.0;                            % Initializing if resetting to different stimuli
        D(1+Nt*(stim-1),:) = 1.0;
        S(1+Nt*(stim-1),:) = 0.0;
%     else
%         r(1+Nt*(stim-1),:) = r(Nt*(stim),:) ;               % Do not initialize if continuing to count stimuli
%         D(1+Nt*(stim-1),:) = D(Nt*(stim),:);
%         S(1+Nt*(stim-1),:) = S(Nt*(stim),:);
    end

    %% Step Through Times
    startInd   = max(trialBin(trial,1),2);
    stopInd    = trialBin(trial,2);
    for i = startInd:stopInd

        I = S(i-1,:)*W+Iapp(i,:) ...        % I depends on feedback (W*S) and applied current
            + sigma*randn(MainStream,1)/sqrt(dt);   % and additional noise

        rinf = rMax./(1.+exp(-(I-Ith)./Iwidth));        % Firing rate curve gives the steady state r
        r(i,:) = rinf + (r(i-1,:)-rinf).*exp(-dt./tauM);  % Update r from the previous timestep

        Dinf = 1./(1.+p0.*r(i-1,:).*tauD);                  % Steady state value of D for Poisson spiking
        D(i,:) = Dinf + ( D(i-1,:)-Dinf).*...
            exp(-dt*(p0.*r(i-1,:)+1./tauD));  % Update with adjusted time constant

        Sinf = sFrac.*p0.*r(i,:).*D(i,:).*tauS./(1.0+sFrac.*p0.*r(i,:).*D(i,:).*tauS); % Steady state value of synaptic gating vatiable assuming vesicle release at a rate p0*r*D
        S(i,:) = Sinf + ( S(i-1,:)-Sinf) .* ...
            exp(-dt*(sFrac.*p0.*r(i,:).*D(i,:)+1./tauS)); % update S with adjusted tau
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
        imagesc(t,1:nCells,r(:,1:end-1)');
        colorbar
        drawnow
    end

%     meanrate = meanrate + r;
%     stdrate = stdrate + r.*r;

end

%% Post-simulation Analysis
