classdef InputStimulus
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
        neurIdentities; % Universal identity vector that labels how many neurons of what types, and in which locations in the vector
        dt;

        % Control of trials
        nTrials=1;

        % Times that stimulus remains on
        trialDuration;
        trialStart;
        trialStop;

        % Hetergeneity, fraction of neurons getting input
        iFrac=1; % Fraction of neurons receiving stimulus
        stimulateInhibitoryCells=false;

        % Number of states stimulus can take, once per trial
        nAmps=1;            % number of states
        opposites=false;    % options, false-none, 'negatives'-negative current, 'population'-partition set halves -- they're all not totally realistic
        opposingPopInh=0; % if provide a fraction, then a fraction of the excitation for the active population is inhibiting opposing population

        % Flag - GPU vector
        useGPU=false;

        % See interface superclass for other defintions
    end

    % METHODS TO GEN STIMULI
    methods
        % ---------------------------------------------------------------
        function this = generateStimuli(this,yokedStim)

            %% Pre-processing steps

            % Establish a random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');

            %% Compute stimulus features
            % Determine the fraction of stimuli that will receive input
            whoStimulate = (rand(r,size(this.neurIdentities))- this.iFrac) ...
                < 0;
            if ~this.stimulateInhibitoryCells
                whoStimulate ...
                    = logical(whoStimulate .* (this.neurIdentities == 0));
            end
            if isequal(this.opposites,'population')
                oppStimulate = find(~whoStimulate);
                oppStimulate = oppStimulate(randperm(numel(oppStimulate)));
                oppStimulate = oppStimulate(...
                    1:...
                    round(this.iFrac*numel(this.neurIdentities))...
                    );
                oppStimulate = ismember(1:numel(this.neurIdentities), oppStimulate);
                oppStimulate ...
                    = logical(oppStimulate .* (this.neurIdentities == 0));
            else
                oppStimulate = zeros(size(this.neurIdentities));
            end

            % Bin up time
            times = 0:this.dt:this.trialDuration*this.nTrials; %#ok<*PROP,*PROPLC>
            this.times = times;

            % Need to calculate start and stops for stimulus
            stimulus = mod(times,this.trialDuration)>=this.trialStart & ...
                mod(times,this.trialDuration)<=this.trialStop;
            % Now replicate into a neuron sized matrix
            stimulus = repmat(stimulus',[1 numel(this.neurIdentities)]);
            % And then filter out neurons unstimulated by the specified
            % fraction
            stimulus(:,~(whoStimulate|oppStimulate)) = 0;

            % acquire the bin numbers of the times at which the stimulus is
            % on. we can do this by looking at the diff'd mod'd times with
            % respect to trial durations. the moment at which the modulus
            % remainder drops.
            diffMod = diff(mod(times,this.trialDuration));
            indStart = find(diffMod < 0); indStart = [1, indStart(1:end-1)];
            indStop = [indStart(2:end)-1 numel(times)];
            assert(numel(indStart) == numel(indStop));

            % Randomly pick stimulus values or compute them from a yoked
            % stimulus.
            % create a random stimulus strength to present the stimulus
            % at!
            if nargin == 1
                trialval = ceil( (this.nAmps)*rand(r,1,this.nTrials) );
            else
                assert(isequal(class(yokedStim),'InputStimulus_Simple'));
                trialval = yokedStim.nAmps-abs(yokedStim.trialval);
            end
            if this.opposites
                stimsign = (round(rand(1,this.nTrials))*2)-1;
                trialval = trialval .* stimsign;
            end
            
            % now for each indStart/indEnd pair, need to randomly assign a
            % strength
            cnt=0; stimulus=cast(stimulus,'double');
            for i = [indStart;indStop]

                cnt=cnt+1;
                
                if isequal(this.opposites,'population')
                    if trialval(cnt) >= 0
                        stimulus(i(1):i(2),whoStimulate) = ... 
                            (stimulus(i(1):i(2),whoStimulate) * trialval(cnt));
                        stimulus(i(1):i(2),oppStimulate) = ... 
                            stimulus(i(1):i(2),oppStimulate) * ...
                            -(this.opposingPopInh * trialval(cnt));
                        stimulus(i(1):i(2),~whoStimulate) = 0;
                    else
                        stimulus(i(1):i(2),oppStimulate) = ...
                            stimulus(i(1):i(2),oppStimulate) * abs(trialval(cnt));
                         stimulus(i(1):i(2),whoStimulate) = ... 
                            stimulus(i(1):i(2),whoStimulate) * ...
                            -(this.opposingPopInh * trialval(cnt));
                        stimulus(i(1):i(2),~oppStimulate) = 0;
                    end
                else
                    stimulus(i(1):i(2),whoStimulate) = ... 
                      (stimulus(i(1):i(2),whoStimulate) * trialval(cnt));
                end

                % store the start and stop indices of a trial
                this.trialBin = [this.trialBin; i(1) i(2)];
            end
            
            % Convert into a scaled stimulus
            stimulus = stimulus * this.I_base;

            %% Output
            this.Iapp = stimulus; % Value of current for each cell
            this.trialval = trialval; % Value of the stimulus for the trial (Useful later for the multiregresion)
            this.stimulus = stimulus(:,find(whoStimulate, 1, 'first')); % Value of stimulus for a cell receiving the input
            if isequal(this.opposites,'population'),
                this.stimulus=this.stimulus-...
                    stimulus(:,find(oppStimulate, 1, 'first')); 
            end

        end
        % ---------------------------------------------------------------
        function [Iapp, trialval, stimulus, trialBin, times] = returnOutputs(this)

            try
                if this.useGPU
                    Iapp = gpuArray(single(this.Iapp));
                    trialval=this.trialval;
                    stimulus=this.stimulus;
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
                trialval=this.trialval;
                stimulus=this.stimulus;
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

        % Stimulus record for a cell receiving the input
        stimulus;

        % Stimulus value on a trial by trial basis, not in time, but by trial
        trialval;
    end


end
