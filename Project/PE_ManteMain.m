% PE_ManteMain -- Parameter Explorer (PE) Script
%
% This script calls a parameter exploration class I wrote that can either
% check all combinations of a particular parameter set or explore all
% possible combinations.
%
% The main benefits :
% (A) Executes each parameter combination in parallel on a separate core,
%       so it explores much faster on a single computer.
% (B) Saves progress, automatically, so that if it stops 12 out of 80
%       simulations deep for any reason, it picks its progress back up,
%       skipping completed sets.
% (C) Contains a framework for automatically placing results into a nested
%       hierarchy of folders that represent the distinct combination per
%       set.
% (D) Records whatever is printed in the console during a sim to respective
%       folder.
% (E) Other constants that currently are not changing parameters can be given,
%       and PE will record them in the project directory, so you never run a sim
%       with a static parameter that you later change, and forget its value.

%% Flags and Shortcut Functions
m2c = @(x) {num2cell(x)};

%% Script and Project Location
% Declare the script we will iterate through, the folder to store resultant
% sims, and (optionally), a session ID to assign to each completed sim (one will
% be randomly created if not).
scriptname      = 'ManteMain';
projectfolder   ='~/Data/Miller/Test';

%% Variables
% Set parameters -- the space over which all value combinations explored
param = struct( ...
       ... Model noise
       'sigma', {{0.1,0.2}}    ...
	);

% Set constants -- this variable space remains the same for each script run
% ...BUT are each actually documented in the containing folder, for the
% purpose that we may later choose to vary them!
consts = struct( ...
     ... Model Resolution/Precision
           'dt',    2e-4, ... Model bin size
           'bin_t', 5e-2, ... How to bin post hoc for analysis
           ... Model noise
           'sigma', 0.1,    ...
       ... ---TASK PARAMETERS---
           ... General Trial Controls
           'nTrials',       7, ...
           'trialDuration', 2,  ... 3, ...
           ... General Stimulus
           'iFrac',         0.15, ...   Randomly select a third of the population for each stimulus to receive
           ... How Opposing Stim Handled (how to compute cells stim fed to)
           'oppSimple',     false, ...
           ... Context Stimulus
           'context_trialStart',    0.1, ...
           'context_trialEnd',      1.4, ...
           ... Color Motion Stimulus
           'color_trialStart',      0.6, ... 0.35, ...
           'color_trialEnd',        1.6,    ... 1.1, ...
           ... Dot Color Stimulus
           'dot_trialStart',        0.6, ... 0.35, ...
           'dot_trialEnd',          1.6, ... 1.1, ...
       ... ---NEURAL NETWORK PARAMETERS---
           ... Number of neurons of each type
           'nInh',      5, ...
           'nExc',      100, ...
           ... Neural Properties
           'rMax0E',    100, ...
           'rMax0I',    200, ...
           'p0E',       0.40, ...
           'p0I',       0.1, ...
           ...Time Constants
           'tausE',     0.025, ...
           'tausI',     0.005, ...
           'tauDbar',   0.2, ...
           'tauDvar',   0.01,     ...
           'tauM',      0.010, ...
           ... Input Properties
           'Imax',      6, ... scales input such that max equals this
           'IthE',      17.75, ... current needed for half maximal input
           'IthI',      20, ...
           'IwidthE',   3, ...
           'IwidthI',   5, ...
           ... Connection Properties
           'WEErecurrent_factor',       200, ...
           'WEEasym_factor',            50,  ...
           'WIE_factor',                -320, ...
           'WEI_factor',                320, ...
           'sigmaWEErec',               0 , ...
           'sigmaWEEasym',              0, ...
           'sigmaIE',                   0, ...
           'pEE',                       0.06 ...
       );
%% Create Explorer, Set Options

% Run the parameter explorer for with reward versus not reward conditions!
Explorer = ParameterExplorer(scriptname,param,...
    'projectfolder',projectfolder,'consts',consts);
% Set options to run with, one specifying that, yes, we will run in
% parallel; the other, that jobs ought be deleted after they're finished
Explorer.useparallel    = true;
Explorer.deletejobs     = true;
Explorer.const2parm     = true;
%Explorer.start = 0.5;
% Set session ID - if not selected, PE will pick one at random
Explorer = Explorer.setSessionID('7');

%% RUN!
Explorer.run();
