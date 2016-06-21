classdef UnitInterface % Abstract Class
    % Template class for defining/organizing code that specifies and readies the
    % parameters descring the collect of neural units in the model
    
    % INPUT PARAMS
   properties
       neurIdentities;
   end
   
   % METHOD -- defines the output parameters, using input parameters, and
   % parameters provided by to the constructor
   methods (Abstract)
       generateOutputParams(this);
   end
   
   % OUTPUT PARAMS
   properties
       % Time Constants
       tauM;
       tauD;
       tauS;
       
       % Synaptic release/binding probs/fractions
       sFrac;
       p0;
       
       % Max Firing rates
       rMax;
   end
end