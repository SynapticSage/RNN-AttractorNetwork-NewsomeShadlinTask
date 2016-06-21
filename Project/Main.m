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
       ... Number of neurons of each type
       'nInh',  5, ...
       'nExc',  40, ...
       ... Model Resolution/Precision
       'dt',    2e-4, ...
       'bin_t', 5e-2, ...
       ... General Trial Controls
       'nTrials',       100, ...
       'trialDuration', 3, ...
       ... Dot Motion Stimulus
       ... Dot Color Stimulus
       ... Context Stimulus
       ... N
       );
end

% Setup neuron counts
neurIdentites=[false(1,params.nExc) true(1,params.nInh)];

%% Initialize Parameter-based Simulation Vectors

% ------------------------------------------------------------------------
% NEURAL VECTORS

% ------------------------------------------------------------------------
% CONECTION VECTORS

% ------------------------------------------------------------------------
% STIMULI VECTORS

%% Execute Simulation

%% Post-simulation Analysis