classdef InputStimulus_Simple
    % This class calculates a the application of a stimulus across trials
    % for a particular stimulus and a set of neurons. It's a simple
    % stimulus in that it's a numbered current with a weighted application
    % to each neuron.
    %
    % TODO : Implement random selection of different stimulus levels!

    % INPUT PARAMS
    properties

        % General stimulus property
        I_base = 1;     % Default scale of the stimulus ..

        % General things that the stimulus has to have acces to to
        % calculate inputs wrt neurons and model time
        neurIdentities; % Universal identity bector that labels how many neurons of what types, and in which locations in the vector
        dt;

        % Control of trials
        nTrials=1;

        % Times that stimulus remains on
        trialDuration;
        trialStart;
        trialStop;

        % Hetergeneity, fraction of neurons getting input
        iFrac=1; % Fraction of neurons receiving stimulus

        % Number of states stimulus can take, once per trial
        nStates=1;      % number of states
        positiveAndNegative=true; % increases the number of states by a factor of two if on

        % Flag - GPU vector
        useGPU=false;

        % See interface superclass for other defintions
    end

    % METHODS TO GEN STIMULI
    methods
        % ---------------------------------------------------------------
        function this = generateStimuli(this)

            %% Pre-processing steps

            % Establish a random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');

            %% Compute stimulus features
            % Determine the fraction of stimuli that will receive input
            whoStimulate = (rand(r,size(this.neurIdentities))- this.iFrac) ...
                < 0;

            % Bin up time
            times = 0:this.dt:this.trialDuration*this.nTrials; %#ok<*PROP,*PROPLC>
            this.times = times;

            % Need to calculate start and stops for stimulus
            stimulus = mod(times,this.trialDuration)>=this.trialStart & ...
                mod(times,this.trialDuration)<=this.trialStop;
            % Convert into a scaled stimulus
            stimulus = stimulus * this.I_base;
            % Now replicate into a neuron sized matrix
            stimulus = repmat(stimulus',[1 numel(this.neurIdentities)]);
            % And then filter out neurons unstimulated by the specified
            % fraction
            stimulus(:,~whoStimulate) = 0;

            % acquire the bin numbers of the times at which the stimulus is
            % on. we can do this by looking at the diff'd mod'd times with
            % respect to trial durations. the moment at which the modulus
            % remainder drops.
            diffMod = diff(mod(times,this.trialDuration));
            indStart = find(diffMod < 0); indStart = [1, indStart];
            indStop = [indStart(2:end)-1 numel(times)];
            assert(numel(indStart) == numel(indStop));

            % now for each indStart/indEnd pair, need to randomly assign a
            % strength
            for i = [indStart;indStop]

                % create a random stimulus strength to present the stimulus
                % at!
                stimStrength = ceil((this.nStates)*rand(r,1,1));
                if this.positiveAndNegative
                    stimSign = round(rand(1))*(-1);
                    if stimSign, stimStrength = stimSign * stimStrength; end; %#ok<BDLGI>
                end
                stimulus(i(1):i(2),whoStimulate) = ...
                  stimulus(i(1):i(2),whoStimulate) * stimStrength;

                % store the start and stop indices of a trial
                this.trialBin = [this.trialBin; i(1) i(2)];
            end

            %% Output
            this.Iapp = stimulus;

        end
        % ---------------------------------------------------------------
        function [Iapp, trialBin, times] = returnOutputs(this)

            try
                if this.useGPU
                    Iapp = gpuArray(single(this.Iapp));
                    trialBin = this.trialBin;
                    times = this.times;
                end
            catch
                this.useGPU=false;
                warning('GPU Memory too small; using CPU');
            end

            if ~this.useGPU
                Iapp = this.Iapp;
                trialBin = this.trialBin;
                times = this.times;
            end
        end
        % ---------------------------------------------------------------
    end

    % OUTPUT PROPERTY
    properties (SetAccess = protected)

        % Specific to this class
        trialBin = [];   % Record of start and stop bin per trial

        % Store a record of experiment time axis
        times;

        % General applied current due to the stimulus to the n units in the
        % simulation.
        Iapp;
    end


end
