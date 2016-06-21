% multi_unit_ISN2
% Test to see bistability in ISN regime with two groups
clear

%% Initialize random number generator
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s);

%% Set up parameters
Nunits = 20;        % No. of cell-groups
Ne = Nunits;        % No. of excitatory rates
Ni = Nunits;        % No. of inhibitory rates

%% Single cell response properties
remax = 50;         % Maximum rate of excitatory cells
rimax = 100;        % Maximum rate of inhibitory cells
rethresh = 20;      % Current for half maximum rate of e-cells
rithresh = 30;      % Current for half maximum rate of i-cells
resigma = 10;       % Range of currents to change e-cell rate
risigma = 10;       % Range of currents to change i-cell rate

% Equations for sigmoidal f-I curves below
fe = @(x) remax./(1+exp(-(x-rethresh)/resigma));    % e-cells
fi = @(x) rimax./(1+exp(-(x-rithresh)/risigma));    % i-cells

%% Connectivity parameters
Wee0 = 1.75;         % E-to-E within unit
Wei0 = 1;         % E-to-I within unit
Wie0 = 0.5;         % I-to-E within unit
Wii0 = 1;           % I-to-I within unit

Weexrand = 0.5/(Nunits-1);    % E-to-E between units (random range)

Weix = 7.5/(Nunits-1);        % E-to-I between units (fixed strength)

%% E-to-I cross-connections are equal strength but probabilistic, with
%  constraint such that each unit provides and receives the same number.
Weiprob = 0.5;                  % Probability of E-to-I connection
Num_ei = floor(Weiprob*(Ne-1)); % Number of E-to-I connections per unit.
Weiconn = zeros(Ne,Ni);         % Weiconn indicates which connections exist

% Use Latin square to equalize connection numbers
a = latinsq_maxdiag(Ne);    % Rows swapped so "Ne" falls along the diagonal
Weiconn = a < Num_ei+1;     % Connectivity matrix is binary

Wee = Wee0*eye(Nunits) + Weexrand*rand(Nunits);
Wei = Weix.*Weiconn + Wei0*eye(Nunits);
Wie = Wie0*eye(Nunits);         % I-to-E connections are local
Wii = Wii0*eye(Nunits);         % I-to-I connections are local

%% Set up simulation paramters
dt = 0.0001;                    % time-step
tmax = 6;                       % maximum time
tvec = 0:dt:tmax;               % vector of time-points
Nt = length(tvec);              % length of time-vector
taue = 0.005;                   % excitatory cell time constant
taui = 0.005;                   % inhibitory cell time constant

noise_I = 150;                  % cell-independent noise added to current

%% Initialize variables
re = zeros(Nt,Ne);              % rates of excitatory cells
ri = zeros(Nt,Ni);              % rates of inhibitory cells
Ie = zeros(Nt,Ne);              % total input current to e-cells
Ii = zeros(Nt,Ni);              % total input current to i-cells

%% Iterate through all time-steps
for i = 2:Nt
    % input currents from e-feedback - i-feedback + noise
    Ie(i,:) = re(i-1,:)*Wee - ri(i-1,:)*Wie  ... 
        + randn(1,Ne)*noise_I*sqrt(dt/taue);
    Ii(i,:) = re(i-1,:)*Wei - ri(i-1,:)*Wii  ...
        + randn(1,Ni)*noise_I*sqrt(dt/taui);
    
    % update firing rates based on f-I curves
    re(i,:) = re(i-1,:)*(1-dt/taue) + dt*fe(Ie(i,:))/taue;
    ri(i,:) = ri(i-1,:)*(1-dt/taui) + dt*fi(Ii(i,:))/taui;
    
end

%% Set up graphics for output
set(0,'DefaultLineLineWidth',3,...
    'DefaultLineMarkerSize',20, ...
    'DefaultAxesLineWidth',3, ...
    'DefaultAxesFontSize',20,...
    'DefaultAxesFontWeight','Bold');

figure()
imagesc(re')        % firing rates of e-cells in heat map
caxis([0 50])
colorbar
xlabel('Time (sec)')
ylabel('Neural Group No.')
c = colorbar;
set(c,'FontSize',20)
ylabel(c, 'Firing Rate (Hz)')
set(gca,'XTick',[20000 40000 60000 80000])
set(gca,'XTickLabel',{'2', '4', '6', '8'})

