function [SingleCondition, PermCondition] =  populationAnalyses(TrialCollection,dt)
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
Xspec = struct('condition',[],'cValue',[],'data',[]);  % vector storing each condition  
Xall = []; 
cnt = 0;
F_single = F(:,:,1); %copy a redundant page from the 3D matrix (each page was a copy for the cell)
for c = 1:size(F_single,1)
  % For each unique value of that condition
  uniques = unique(F_single(c,:));
  for u = uniques
    cnt = cnt+1;
    
    % Identify each trial bearing the value
    trials = u == F_single(c,:);
    Xspec(cnt).condition = c;
    Xspec(cnt).cValue = u;
    popResp = populationAverage(r,trials);
    Xspec(cnt).data = popResp;
    Xall = [Xall; popResp];
  end
end
% Compute the number of bins in a condition
conditionTLength = size(popResp,1);

%% Permutation population analysis  - permutation stimulus averages, 
% e.g. context 1 then motion 2, being a permutation
%  ... commented out for now! need to create a separate beta with an
%  element for each of the permutations

Yspec = struct('condition',[],'cValue',[],'data',[]);  % vector storing each condition  
Yall = []; 
conditions = unique(F_single','rows');
nConditions = size(conditions,1);
permPopResp = zeros(nTime,nCells,nConditions);
uPerm = unique(F_single','rows');
cnt = 0;
for u = uPerm'
    cnt=cnt+1;
    % find matching trials!
    [~,trials] = ismember(u',F_single','rows');
    % Setup store trial characteristics
    Yspec(cnt).cValue = u';
    popResp = populationAverage(r,trials);
    Yspec(cnt).data = popResp;
    Yall = [Yall; popResp];
end

%% Acquire Targeted Dimensionality Reduction Matrices

% Obtain denoising matrices
    % First for condition-specifc population averages
    xD = denoisingMatrix(Xall);
    % Second for condition permutation population averages
     yD = denoisingMatrix(Yall);
    
% Obtain regression subspace
    % First for condition-specifc population averages
    orthX = regressionSubspace(r,beta,xD);
%     % Second for condition permutation population averages
     orthY = regressionSubspace(r,beta,yD);

%%  Project 

% Obtain projection onto the three axes
for c = 1:numel(Yspec)
    Yspec(c).proj = Yspec(c).data*orthY;
end

% Obtain projection onto the three axes
for c = 1:numel(Xspec)
    Xspec(c).proj = Xspec(c).data*orthX;
end
     
%% Package up all the statistical measures into struct
%

SingleCondition.specifics = Xspec;
SingleCondition.all = Xall;
SingleCondition.orthogConditionAxes = orthX;
SingleCondition.denoise = xD;

PermCondition.specifics = Yspec;
PermCondition.all = Yall;
PermCondition.orthogConditionAxes = orthY;
PermCondition.denoise = yD;


%% 
    function [orthAxes, Bpac_max] = regressionSubspace(r,beta,D)
        % Output, orthogonal axes Nunits X NConditions, and also a set of
        % non-orthogonalized "axes", because is informative of overlap in
        % population code.
        
        % First, we submit beta to the same denoising as can be done to a
        % populatino vector, across its time and cellular dimensions.
        % ... beta needs to be permuted so that it can be mutiplied
        beta = permute(beta, [3 2 1]);
        Bpca = zeros( size(beta) );
        for i = 1:size(beta,3) % over number of condtions
            Bpca(:,:,i) = D * beta(:,:,i);
        end 
        
        % Now we determine the time at which these beta vectors have
        % maximal norm! We can then derive time-independent, de-noised
        % regression vectors.
         Bpca_norms = sqrt(sum(Bpca.^2,1));
        [~,t_max_v] = max(Bpca_norms,[],2);
        t_max_v = squeeze(t_max_v);
        
        Bpca_max = zeros( size(Bpca,1), size(Bpca,3) );
        for i = 1:size(t_max_v,3)
            Bpca_max(:,i) = Bpca(:,t_max_v(i),i);
        end
        
        
        % Now we need to orthogonalize them with respect to eachother, to
        % obtain the orthogonal components we're going to use to describ
        % the population
        [Q,R] = qr(Bpca_max);
        
        % Slice out the orthogonal components
        orthAxes = Q(:,1:size(beta,3));
        
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
