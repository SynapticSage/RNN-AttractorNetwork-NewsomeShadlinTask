function populationAnalyses(TrialCollection,dt)
% Author: Ryan Y.
% Function executes much of the analysis done in the Mante et al. 2015
% paper. In order for this to work well, run a huge number of trials, if
% posisble, such that trial conditions are repeated a number of times!

r = TrialCollection.r;
F = TrialCollection.F;

nCells = size(r,3);
nTime = size(r,2);
nTrial = size(r,1);

%% GLM

% To perform GLM, we will take inverse(F*F')F*r, on each matrix indiced in the
% third dimension across neurons.

% If mmx-libraries exist on this machine, initiate the ultra fast version of
% this computation, because these are ixjxnCell, and mmx allows us to vectorize
% matrix multiplication
usemmx = exist('mmx.m','file');
try
if usemmx

  squareF = mmx('square',F,[]); % performs A*A' on each dim
  Fr = mmx('mult',F,r); % performs F*r on each dim
  beta = mmx('backslash',squareF,Fr); % performs (squareF)\Fr per dim
end
catch
% If an exception is caught, user probably hasn't built the mex code for
% their system, so we should fall back on vanilla matlab commands.
usemmx = false;
end
% Otherwise, we have to do this the slow!, "complicated" way.
if ~usemmx
  beta = zeros(size(F,1),size(r,2),nCells);
  squareF = zeros(size(F,1),size(F,1),nCells);
  Fr = zeros(size(F,1),size(r,2),nCells);
  for c = 1:nCells
    squareF(:,:,c) = F(:,:,c)*F(:,:,c)';
    Fr(:,:,c) = F(:,:,c)*r(:,:,c);
    beta(:,:,c) = squareF(:,:,c)\Fr(:,:,c);
  end
end

% What we now have is a beta that predicts the cells response based on the task
% characteristics at each time in the task. Beta tells us how each variable
% affects that cell per time point.

%% Population analysis - single stimulus averages
nConditions = size(F,1);
X_spec = struct('condition',[],'cValue',[],'data',[]);  % vector storing each condition  
X_all = []; 
cnt = 0;
F_single = F(:,:,1); %copy a redundant page from the 3D matrix (each page was a copy for the cell)
for c = 1:size(F_single,1)
  % For each unique value of that condition
  uniques = unique(F_single(c,:));
  for u = uniques
    cnt = cnt+1;
    
    % Identify each trial bearing the value
    trials = u == F_single(c,:);
    X_spec(cnt).condition = c;
    X_spec(cnt).cValue = u;
    popResp = populationAverage(r,trials);
    X_spec(cnt).data = popResp;
    X_all = [X_all; popResp];
  end
end
% Compute the number of bins in a condition
conditionTLength = size(popResp,1);

%% Population analysis  - permutation stimulus averages
Y_spec = struct('condition',[],'cValue',[],'data',[]);  % vector storing each condition  
Y_all = []; 
conditions = unique(F_single','rows');
nConditions = size(conditions,1);
permPopResp = zeros(nTime,nCells,nConditions);
uPerm = unique(F_single','rows');
for u = uPerm'
    % find matching trials!
    [~,trials] = ismember(u',F_single','rows');
    % Setup store trial characteristics
    Y_spec(cnt).cValue = u';
    popResp = populationAverage(r,trials);
    Y_spec(cnt).data = popResp;
    Y_all = [Y_all; popResp];
end

%% Targeted Dimensionality Reduction

% Obtain denoising matrices
    % First for condition-specifc population averages
    xD = denoisingMatrix(X_all);
    % Second for condition permutation population averages
    yD = denoisingMatrix(Y_all);
    
% Obtain regression subspace
    % First for condition-specifc population averages
    xBpca = regressionSubspace();
    % Second for condition permutation population averages
    yBpca = regressionSubspace();


%% 
    function regressionSubspace(r,beta,D)
        % First, we submit beta to the same denoising as can be done to a
        % populatino vector, across its time and cellular dimensions.
        
        % Now we determine the time at which these beta vectors have
        % maximal norm! We can then derive time-independent, de-noised
        % regression vectors.
        
        % Now we need to orthogonalize them with respect to eachother, to
        % obtain the orthogonal components we're going to use to describ
        % the population
        
        % Slice out the orthogonal components
        
        % These are the vectors we will return, such that we can the
        % project population averages per condition onto them!
        
    end

    function X = populationAverage(r,trials)
        % Shortcut function for population averaging. It's also convenient
        % because a change here (e.g. addition of smoothing) will affect
        % both sections that do something like this above.
        
        X = mean(r(trials,:,:),1);
        X = squeeze(X);
        
%         plot(X);
        
        % % UNCOMMENT IF USE POISSON MODEL INSTEAD OF AN ALREADY SMOOTH FIRING RATE
        % % Create guassian kernel of of 40ms
        % smoothWidth = 40e-3;
        % kt = [-4:dt:4]';
        % kernel = exp(-((kt).^2/(2*smoothWidth)^2));
        % % Permute time axis into first dim if not already
        % % Replicate kernel over neuron number and condition numbers
        % kernel = repmat(kernel,1,nCells,nConditions);
        % Xsmooth = conv(X,kernel,'same');
    end

    function [D] = denoisingMatrix(popMat,varargin)
        % This subfunction carries out the Targeted Dimensionality
        % Reduction technique in the Mante et al. 2015 paper.
        %
        % Outputs ... denoising matrix
        
        nComponents = 12;  % Default number of components to accept into the denoising vector.
        for v = 1:2:numel(varargin)
            switch varargin{v}
                % nComponents can be provided to change the default number
                % of smoothing components.
                case 'nComponents', nComponents = varargin{v+1};
            end
        end
       
        % Obtain prominent pca components
        [comp,score] = pca(popMat);
        comp =  comp(:,1:nComponents); 
        score = score(:,1:nComponents);
        % Construct the Denoising matrix, D
        D = zeros(nCells,nCells);
        for c = comp
            D = D + c*c';
        end
        
    end

end
