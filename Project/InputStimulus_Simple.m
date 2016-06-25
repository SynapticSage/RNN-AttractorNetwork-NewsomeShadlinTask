classdef InputStimulus_Simple < StimuliInterface
    % This class calculates a the application of a stimulus across trials
    % for a particular stimulus and a set of neurons. It's a simple
    % stimulus in that it's a numbered current with a weighted application
    % to each neuron. 
    %
    % TODO : Implement random selection of different stimulus levels!
    
    % INPUT PARAMS
    properties 
        % Strength Levels
        
        % Control of trials
        nTrials=1;
        
        % Times that stimulus remains on
        trialDuration;
        trialStart;
        trialStop;
        
        % Hetergeneity, fraction of neurons getting input
        iFrac=1; % Fraction of neurons receiving stimulus
        
        % Number of states stimulus can take, once per trial
        nStates=1;
        
        % Flag - GPU vector
        useGPU=true;
        
        % See interface superclass for other defintions
        
    end
    
    % METHODS TO GEN STIMULI
    methods
        % ---------------------------------------------------------------
        function this = generateStimuli(this,dt)
            
            %% Pre-processing steps
            % Create an alias for this object
            t=this;
            
            % Establish a random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');
            
            %% Compute stimulus features
            % Determine the fraction of stimuli that will receive input
            whoStimulate = (rand(r,size(t.neurIdentities))- t.iFrac) < 0;
            
            % Bin up time
            times = 0:t.dt:t.trialDuration*t.nTrials;
            
            % Need to calculate start and stops for stimulus
            stimulus = mod(times,t.trialDuration)>=t.trialStart & ...
                mod(times,t.trialDuration)<=t.trialStop;
            % Convert into a scaled stimulus
            stimulus = stimulus * t.I_base;
            % Now replicate into a neuron sized matrix
            stimulus = repmat(stimulus,[numel(this.neurIdentities) 1]);
            % And then filter out neurons unstimulated by the specified
            % fraction
            stimulus(:,~whoStimulate) = 0;
            
            % acquire the bin numbers of the times at which the stimulus is
            % on. we can do this by looking at the diff'd mod'd times with
            % respect to trial durations. the moment at which the modulus
            % remainder drops.
            diffMod = diff(mod(times,t.trialDuration));
            indStart = find(diffMod < 0) + 1;
            indEnd = find(diffMod < 0);
            assert(numel(indStart) == numel(indEnd));
            
            
            %% Output
            t.Iapp = stimulus;
            this=t;
            
        end
        % ---------------------------------------------------------------
        function [Iapp, trialBin] = returnOutputs(this)
            if this.useGPU
                Iapp = gpuArray(this.Iapp); 
                trialBin = gpuArray(this.trialBin);
            else
                Iapp = this.Iapp; 
                trialBin = this.trialBin;
            end
        end
        % ---------------------------------------------------------------
    end
    
    % OUTPUT PROPERTY
    properties (SetAccess = protected)
        % See interface for other definitions
        
        % Specific to this class
        trialBin;   % Record of start and stop bin per trial
    end

    
end