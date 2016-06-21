% many_tuned_count_mems_50h6.m
% Rate-model code for a network with many stable states.

%% Pre-processing
clear
%close all
clf

%% Instantiate Properties of Simulation
% DIFFERENT RANDOMNESS STREAMS to ensure different aspects not correlated
% 2 random streams, with s2 used for cell and network structure while the
% other, s1, is used for temporal noise
[s1, s2 ] = ...
    RandStream.create('mrg32k3a','NumStreams',2,'seed',sum(100*clock));
% [s1 s2 ] = RandStream.create('mrg32k3a','NumStreams',2,'seed',2);

% --- SIMULATION: PRECISION ---
dt      = 0.0001;
bint    = 0.050;          % what does this control?
nbin    = round(bint/dt)  % number of bins in a time bin? or total?

% How long to run and how many bins per second
tmax    = 1.5;            % seconds
Nsec1   = round(0.75*tmax/dt); % 3/4 of the total time, split into dt sized bins

% --- TRIALS: REPITITION ---
Ntrials         = 10;   % total number of trials
trial_reset     = 0;    % whether to reset after trial?
dtime           = true;    % is this a flag?
Num_to_average  = 10;

% --- ENTIRE TRIAL: TIME ---
t   = 0:dt:tmax*(Ntrials+1);      % Time vector
Nt  = round(tmax/dt);            % Total number of bins in the complete time

% --- CELLS ---
Ncells      = 51;         % Number of neurons in code

rmax0E      = 100;        % Maximum firing rate
rmax0I      = 200;        % Maximum firing rate
rmax        = [rmax0E*ones(1,Ncells-1) rmax0I]; % vector of rmaxes for all neurons
Ith0        = 6.4;             % Current needed for half-maximal firing
IsigmaTh    = 0;           % Range of thresholds across all cells
Iwidth0     = 1;           % How far from threshold you need to be for rate to change

Wrecurrent  = 85;          % Strength of connection from one group to itself
sigmaWEE    = 0;           % Sigma of the excitatory inputs, in this case 0, indicating a linear response
Wasym       = 20/(Ncells-1);          % Strength of connection from one group to another
sigmaWasym  = 2*Wasym;
WIE         = -300;
sigmaIE     = 0;
WEI         = 250/(Ncells-1);

tausE   = 0.050;
tausI   = 0.005;
taus    = [tausE*ones(1,Ncells-1) tausI];       % Time constant for synapse
taud    = 0.500;    % Time constant for depression
taum    = 0.010;    % Time constant for change in firing rate
p0E     = 1;        % Base release probability of vesicles
p0I     = 0.1;
p0      = [p0E*ones(1,Ncells-1) p0I];
sfrac   = 1;        % Maximum proportion of receptors activated by a vesicle release

dI          = 0.2;
stimtime    = 0.25; % the constant stim time
dstimtime   = 0.005;% the actual different in stimulation time

%% Calculate Stimulus/Stimulation Properties
% Specifies whether or not to simply concatonate the current step onto the
% previous current step, or to simply apply the same stimulus each time
if trial_reset && ~dtime
    Iapp0 = dI:dI:dI*Ntrials;                                               % Applied current steps
else
    Iapp0 =  1*ones(1,Ntrials);
end

if dtime
    trialstimtime = dstimtime   * (1:Ntrials);
else
    trialstimtime = stimtime    * ones(1,Ntrials);
end

Iappsigmafrac = 0.0;                                                        % degree to which applied currents are different
IappHetero = 1 + Iappsigmafrac*(rand(s2,1,Ncells)-0.5);
Istart = 0.1;                                                               % Time to start applied current

sigma = 0.002;                                                              % standard deviation of noise in current

Iapp = zeros(length(t),Ncells);                                             % Applied current for each cell at all time points

for trial = 1:Ntrials
    Iend = Istart + trialstimtime(trial);                                   % Time to finish appplied current
    imin = min(round(Istart/dt)+1,length(t));                               % Index of time for applied current
    imax = min(round(Iend/dt)+1,length(t));                                 % Index of time to end applied current

    for b = imin:imax
        Iapp(b+Nt*(trial),:) = Iapp0(trial)*IappHetero;
    end
end
Iapp(:,Ncells) = 0.0;

%% Current Threshold And Width

Ith     = Ith0      * ones(1,Ncells) + IsigmaTh*(rand(s2,1,Ncells)-0.5);
Iwidth  = Iwidth0   * ones(1,Ncells);

Ith(Ncells)     = 12.0;                                                     % Readout cell has a different threshold
Iwidth(Ncells)  = 3.0;                                                      % Readout cell has a different steepness of slope

%% Setting up cell connections

W = zeros(Ncells,Ncells);                                                   % W is the weight matrix: strength of connections
for cell1 = 1:Ncells
    W(cell1,:)      = Wasym      + sigmaWasym*(rand(s2,1,Ncells)-0.5);      % Set all connections from cell1 to cell2 to Wasym
    W(cell1,cell1)  = Wrecurrent + sigmaWEE*(rand(s2,1)-0.5);               % This is the strength of self-connections and replaces the other value
end
W(:,Ncells) = WEI;                                                          % Connection strength to I-cell from other cells
W(Ncells,:) = WIE + sigmaIE*(rand(s2,1,Ncells)-0.5);
W(Ncells,Ncells) = 0.0;                                                     % I-to-I strength

%% Initialize mean and standard deviation trackers

meanrate        = zeros(length(t),Ncells);
stdrate         = zeros(length(t),Ncells);
mresponse1      = zeros(Ntrials,Ncells,Num_to_average);
meanresponse1   = zeros(Ntrials,Ncells);
sdresponse1     = zeros(Ntrials,Ncells);

%% Trial Computations
tic
for averagenum = 1:2*Num_to_average
    fprintf('Average Number %d\n',averagenum);

    r = zeros(length(t),Ncells);    % Firing rate for each cell at all time points
    D = zeros(length(t),Ncells);    % Depression variable for each cell at all time points
    S = zeros(length(t),Ncells);    % Synaptic gating variable for each cell at all time points

    % Simulate trial
    for trial = 0:Ntrials

        % What we should intitialize the first time point of the tria to...
        if  trial_reset || trial == 0
            r(1+Nt*(trial),:) = 0.0;
            D(1+Nt*(trial),:) = 1.0;
            S(1+Nt*(trial),:) = 0.0;
        else
            r(1+Nt*(trial),:) = r(Nt*(trial),:) ;
            D(1+Nt*(trial),:) = D(Nt*(trial),:);
            S(1+Nt*(trial),:) = S(Nt*(trial),:);
        end

        %% Step Through Time
        % Now calculate the value at each time bin in the trial
        for b = 2+Nt*(trial):Nt*(trial+1)       % Now integrate through time
            
            % Calculate the currents all cells will be receiving
            I = S(b-1,:)*W+Iapp(b,:) ...        % I depends on feedback (W*S) and applied current
                + sigma*randn(1)/sqrt(dt);      % and additional noise
            
            % S(:,i-1) is the vector of synaptic gating from the previous
            % time step for all cells This gets multiplied by the weight
            % matrix to give total feedback current.

            % Firing rate curve gives the steady state r
            rinf = rmax./(1.+exp(-(I-Ith)./Iwidth));
            % Update r from the previous timestep
            r(b,:) = rinf + (r(b-1,:)-rinf)*exp(-dt/taum);

            % Steady state value of D for Poisson spiking
            Dinf = 1./(1.+p0.*r(b-1,:)*taud);
            % Update with adjusted time constant
            D(b,:) = Dinf + (D(b-1,:)-Dinf).* ...
                exp(-dt*(p0.*r(b-1,:)+1/taud));

            % Steady state value of synaptic gating vatiable assuming vesicle release at a rate p0*r*D
            Sinf = sfrac*p0.*r(b,:).*D(b,:).*...
                taus./(1.0+sfrac*p0.*r(b,:).*D(b,:).*taus);
            % update S with adjusted tau
            S(b,:) = Sinf + (S(b-1,:)-Sinf).* ...
                exp(-dt*(sfrac*p0.*r(b,:).*D(b,:)+1./taus));
        end

        %% Compute trial statistics
        if trial > 0
            if averagenum <= Num_to_average
                meanresponse1(trial,:) = meanresponse1(trial,:) + ...
                    mean(r(b-Nsec1:b,:));
                sdresponse1(trial,:) = sdresponse1(trial,:) + ...
                    mean(r(b-Nsec1:b,:)).*mean(r(b-Nsec1:b,:));
            else
                mresponse1(trial,:,averagenum-Num_to_average) =  ...
                    mean(r(b-Nsec1:b,:));
            end
        end

    end
    
    figure(1)
    imagesc(r(:,1:end-1)')
    colorbar
    drawnow

    meanrate = meanrate + r;
    stdrate = stdrate + r.*r;

end
toc;

%% Post-simulation Analysis: Mean Responses

meanresponse1 = meanresponse1/Num_to_average;
sdresponse1 = ...
    sqrt(sdresponse1/Num_to_average-meanresponse1.*meanresponse1);

shiftmresponse1 = meanresponse1 - ones(Ntrials,1)*mean(meanresponse1);
figure(1);
imagesc(meanresponse1');

%% Post-simulation Analysis: Correlations

[scorrr1r1, spvalr1r1]  = corr(shiftmresponse1',shiftmresponse1')
[corrr1r1, pvalr1r1]    = corr(meanresponse1',meanresponse1')

figure(2)
subplot(1,2,1)
imagesc(scorrr1r1)
subplot(1,2,2)
imagesc(corrr1r1)

%% Post-simulation Analysis: Confusibility

bintvec = 0:bint:tmax*(Ntrials+1)-bint;
binmeanrate = zeros(length(bintvec),Ncells);
binstdrate = zeros(length(bintvec),Ncells);

for b = 1:length(bintvec)
    binmeanrate(b,:) = mean(meanrate((b-1)*nbin+1:b*nbin,:));
    binstdrate(b,:) = mean(stdrate((b-1)*nbin+1:b*nbin,:));
end


binmeanrate = binmeanrate/Num_to_average;
binstdrate = binstdrate/Num_to_average - binmeanrate.*binmeanrate;

maxcorrs = zeros(Ntrials,Num_to_average);
for b = 1:Num_to_average
    trialcorr = corr(squeeze(mresponse1(:,:,b))',meanresponse1');
    [dummy, maxcorrs(:,b)] = max(trialcorr);
end

confusibility = hist(maxcorrs',1:10);
confusibility = confusibility/Num_to_average;
figure(3)
imagesc(confusibility)

clear r D S I t meanrate stdrate

save many_count_mems_50_h6_extra
