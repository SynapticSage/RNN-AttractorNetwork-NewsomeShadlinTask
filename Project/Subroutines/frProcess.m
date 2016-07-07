function [firingrate] = frProcess(firingrate,reqSize,varargin)
  % Processes trial information to compute measures like in the Mante, Susillo, Abott paper
  %
  % INPUT:
  %         firingrate -- firing rate for a single trial, times x cells
  %         trialconditions -- vector of trial conditions
  %
  %
  %% Optional Section
  gpuEnable=false;
  for v = 1:2:numel(varargin)
    switch varargin{v}
    case 'gpuEnable', gpuEnable = varargin{v+1};
    end
  end

%% Process
  % First need to create a firing rate that's normalized for that particular time and trial
  meanr = mean(firingrate,1);
  stdr = std(firingrate,1);
  firingrate = (firingrate - meanr)./stdr; temp=firingrate;
  % Need to reshape so that cell count points into 3rd dimension
  firingrate = permute(firingrate,[3 1 2]);
  
%% Ensure size matches required size
% Sometimes its off by a single bin, other times not
if size(firingrate,2) < reqSize(2)
    difference = reqSize(2)-size(firingrate,2);
    firingrate(:,end+1:end+difference,:) = firingrate(:,end,:);
    assert(difference<3,'Bin size mismatch! Check code.');
elseif size(firingrate,2) > reqSize(2)
    difference = size(firingrate,2)-reqSize(2);
    firingrate(:,end-difference+1:end,:) = [];
    assert(difference<3,'Bin size mismatch! Check code.')
end
 

  %% Send to GPU?
  if gpuEnable
    firingrate = gpuArray(firingrate);
  end

end
