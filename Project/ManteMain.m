% Main for the script that implements an attractor based neural firing
% model on a stimulus regime that's supposed to mimic the Mante task.

%% Pre-processing and Basic Definitions

clear; close all;

% Whether to use GPU
useGPU = false;
if useGPU
    reset(gpuDevice(1)); % Resets gpu memory, if anything is in it.
end

% shortcut lambda functions
normalize   = @(x) x./mean(x);
zscore      = @(x) mean(x)./std(x);
dsamp       = @(x) downsample(x,10);

%% General Flags

% ParameterExplorer controlled mode
PE_mode = false;
% Whether to reset network after each trial
trialReset = false;
% Turn on additional plotting
figures.on = true;
figures.showStimuli = true;

%% General Paramters
% The lines below this comment provide an optional entry point for my
% ParameterExplorer class.
if ~(exist('params','var') && PE_mode)
   params = struct( ...
           ... Model Resolution/Precision
           'dt',    2e-4, ... Model bin size
           'bin_t', 5e-2, ... How to bin post hoc for analysis
           ... Model noise
           'sigma', 0.1,    ...
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
           'color_trialEnd',        1.1, ...
           ... Dot Color Stimulus
           'dot_trialStart',        0.35, ...
           'dot_trialEnd',          1.1, ...
       ... ---NEURAL NETWORK PARAMETERS---
           ... Number of neurons of each type
           'nInh',      25, ...
           'nExc',      50, ...
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
           'WEErecurrent_factor',      210, ...
           'WEEasym_factor',    35,  ...
           'WIE_factor',        -320, ...
           'WEI_factor',        320, ...
           'sigmaWEErec',       0, ...
           'sigmaWEEasym',      0, ...
           'sigmaIE',           0, ...
           'pEE',               0.0350 ...
       );
   savedir = '~/Data/Miller/Untitled';
else
    % params and projectfolder already exist provide by PE, now derive
    % savedir from set
    savedir = ParameterExplorer.savelocation(params,...
        'projectfolder',projectfolder);
end

fprintf('SaveLocation: %s\n',savedir);

% Setup neuron counts
neurIdentites=[false(1,params.nExc) true(1,params.nInh)];
nCells = params.nExc + params.nInh;

% Setup random stream for main script (different than random streams that
% are operating to construct stimuli, connections, and neuron properties,
% respectively)
MainStream = RandStream.create('mrg32k3a','NumStreams',1,'seed','shuffle');

% Assign param variables to main space variables that will be directly
% called in the simulation
dt = params.dt;
sigma = params.dt;

%% Initialize Parameter-based Simulation Vectors
% ------------------------------------------------------------------------
% NEURAL VECTORS
NP = NeuronProperties; % Encapsulated code for setting up overall neuron/synaptic vectors
    % First, we set the properties to build the the output properties
    NP.neurIdentities = neurIdentites;
    NP.rMax0E=params.rMax0E;        NP.rMax0I=params.rMax0I;
    NP.p0E=params.p0E;              NP.p0I=params.p0I;
    NP.tausE=params.tausE;          NP.tausI=params.tausI;
    NP.tauDbar=params.tauDbar;      NP.tauDvar=params.tauDvar;
    NP.tauMbar=params.tauM;
    NP.IwidthI =  params.IwidthI;   NP.IwidthE = params.IwidthE;
    NP.Ith0E =  params.IthE;        NP.Ith0I = params.IthI;

    % Then, we invoke the generation method to create them
    NP = NP.generateOutputParams;

    % last return the gpu vectors for the simulation
    [tauM,tauD,tauS,p0,sFrac,rMax,Iwidth,Ith] = NP.returnOutputs();
    clear NP;

% ------------------------------------------------------------------------
% CONECTION VECTORS

% Setup W matrix constructor
Connect = ConnectionProperties; % Creates W vector
    % First, we set the properties to build the W matrix
    Connect.neurIdentities  = neurIdentites;
    Connect.Rec.EE          = params.WEErecurrent_factor;
    Connect.Base.EE         = params.WEEasym_factor;
    Connect.Base.IE         = params.WIE_factor;
    Connect.Base.EI         = params.WEI_factor;
    Connect.Sigma.Rec.EE    = params.sigmaWEErec;
    Connect.Sigma.EE        = params.sigmaWEEasym;
    Connect.Sigma.IE        = params.sigmaIE;
    Connect.P.EE            = params.pEE;

    % then, we invoke the generation method to create W
    Connect = Connect.generateConnections();
    Connect.visualize(savedir);

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
        [Iapp.context, trials.context, trials.raw.context, trials.bin, t] = ...
          Context.returnOutputs();
    % (2) Setup Dot Stimulus
    Dots = GeneralStim;
        Dots.trialStart     = params.dot_trialStart;
        Dots.trialStop      = params.dot_trialEnd;
        Dots.nStates        = 3;
        Dots=Dots.generateStimuli();
        [Iapp.dots, trials.dots, trials.raw.dots] = Dots.returnOutputs();
    % (3) Setup Color Stimulus
    Color = GeneralStim;
        Color.trialStart    = params.color_trialStart;
        Color.trialStop     = params.color_trialEnd;
        Color.nStates       = 3;
        Color=Color.generateStimuli();
        [Iapp.color, trials.color, trials.raw.color] = Color.returnOutputs();

    % Clear uneeded variables from workspace
    clear Context Dots Color GeneralStim;

   % Racast Iapp as the superposition of the inputs influence on cells
   Iapp = Iapp.context + Iapp.color + Iapp.dots;

%% Execute Simulation

r = zeros(length(t),nCells);    % Firing rate for each cell at all time points
D = zeros(length(t),nCells);    % Depression variable for each cell at all time points
S = zeros(length(t),nCells);    % Synaptic gating variable for each cell at all time points
if useGPU
        r=gpuArray(r);
        D=gpuArray(D);
        S=gpuArray(S);
end
for tr = 1:params.nTrials

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
    startInd   = max(trials.bin(tr,1),2);
    stopInd    = trials.bin(tr,2);
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

    %% Calculate post-trial measures
    
    

    %% Post-trial plotting
    if figures.on
        
        figure(1)
        figures.nSubplot = 2;
        
        subplot(figures.nSubplot,1,figures.nSubplot-1);
        imagesc(t,1:nCells-params.nInh,r(:,1:end-params.nInh)'); axis tight;
        colorbar; title('Exc');
        subplot(figures.nSubplot,1,figures.nSubplot);
        imagesc(t,nCells-params.nInh+1:nCells,r(:,end-params.nInh+1:end)'); axis tight;
        colorbar; title('Inh');
        
        drawnow;
    end
end

%% Post-trial plotting
if figures.on
    
    figure(1)
    figures.nSubplot = 2;

    if figures.showStimuli
        figures.nSubplot=figures.nSubplot+1;
        subplot(figures.nSubplot,1,1); hold on;
        
        p1=plot(dsamp(t),dsamp(trials.raw.context),'linewidth',2);
        p2=plot(dsamp(t),dsamp(trials.raw.dots),'linewidth',2);
        p3=plot(dsamp(t),dsamp(trials.raw.color),'linewidth',2);
        
        axis tight;
        
        plot_digital(dsamp(t),dsamp(trials.raw.context),...
            'color',get(p1,'color'));
        plot_digital(dsamp(t),dsamp(trials.raw.dots),...
            'color',get(p2,'color'));
        plot_digital(dsamp(t),dsamp(trials.raw.color),...
            'color',get(p3,'color'));
        
        legend([p1 p2 p3],'Context','Dots','Color');
        
    end

    subplot(figures.nSubplot,1,figures.nSubplot-1);
    imagesc(t,1:nCells-params.nInh,r(:,1:end-params.nInh)'); axis tight;
    colorbar; title('Exc');
    subplot(figures.nSubplot,1,figures.nSubplot);
    imagesc(t,nCells-params.nInh+1:nCells,r(:,end-params.nInh+1:end)'); axis tight;
    colorbar; title('Inh');
    
end

%% Post-simulation statistics
meanrate  = mean(r,1);
stdrate   = std(r,1);

fprintf('\n--------\nMean R:\t');   fprintf('%9.3f ',meanrate);
fprintf('\n--------\nStd R:\t');    fprintf('%9.3f ',stdrate);
fprintf('\n--------\n');

%% Post-simulation Analysis: Multiple Linear Regress

% Construct unit matrices


% Regress trial by trial, motion, color, and context

%

%% Post-simulation Analysis:
