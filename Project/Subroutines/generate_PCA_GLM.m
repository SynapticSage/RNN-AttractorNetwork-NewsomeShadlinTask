function generate_PCA_GLM(tt,dt)

nCells = size(r,3);
nTime = size(r,2);
nTrial = size(r,1);

%% GLM
r = tt.r;
F = tt.trialMat;

% To perform GLM, we will take inverse(F*F')F*r, on each matrix indiced in the
% third dimension across neurons.

% If mmx-libraries exist on this machine, initiate the ultra fast version of
% this computation, because these are ixjxnCell, and mmx allows us to vectorize
% matrix multiplication
if exist('mmx.m','file')

  squareF = mmx('square',F,[]); % performs A*A' on each dim
  Fr = mmx('mult',F,r); % performs F*r on each dim
  beta = mmx('backslash',squareF,Fr); % performs (squareF)\Fr per dim

% Otherwise, we have to do this the slow and "complicated" way.
else

  beta = zeros(size(F,1),size(r,2),nCells);
  squareF = zeros(size(F,1),size(r,2),nCells);
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

%% Population analysis
popResp = zeros(nTime,nCells,nConditions);
F_single = F(:,:,1); %copy the redundant information from the matrix
for C = F_single
  % For each unique value of that condition
  C_unique = unque(C);
  for c = C_unique
    % Generate the population average over-time
  end
end

% % UNCOMMENT IF USE POISSON MODEL INSTEAD OF AN ALREADY SMOOTH FIRING RATE
% % Create guassian kernel of of 40ms
% smoothWidth = 40e-3;
% kt = [-4:dt:4]';
% kernel = exp(-((kt).^2/(2*smoothWidth)^2));
% % Permute time axis into first dim if not already
% % Replicate kernel over neuron number and condition numbers
% kernel = repmat(kernel,1,nCells,nConditions);
% smoothedPopResp = conv(popResp,kernel,'same');

%%

end
