function [firingrate] = trialprocess(firingrate,varargin)
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
  firingrate = (firingrate - meanr)./stdr;

  % Need to output a 1D representation of the firing rate to match the form in the paper
  firingrate = reshape(firingrate,[],1);

  %% Send to GPU?
  if gpuEnable
    firingrate = gpuArray(firingrate);
  end

end
