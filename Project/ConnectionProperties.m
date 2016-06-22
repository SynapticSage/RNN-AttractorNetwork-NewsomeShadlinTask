classdef ConnectionProperties < ConnectionInterface
    % Implementation of connections class, that determines weights for our
    % neural units. A different implementation is called for 2D topography
    % of cells.
    %
    % Honestly, there could an even more flexible implementation, that
    % doesn't presume which connection types have randomness, which don't,
    % which have asmmetric versus recurrent, which don't, and so on. Nor
    % assume there are just two cell populations, namely,
    % excitatory/inhibitory. Which I could see being a huge time saver for
    % trying out other model styles, rapidly.
    
    %% INPUT PROPERTIES
    % Specific to this implementation of the the connections interface
    % class.
    properties (Access = public)  
        % Parameters that specify the asymmetric and reccurent connection
        % properties, and as well as standard ...
        Wee_asym; 
        Wee_recurrent; 
        Wie_standard;
        Wei_standard;
        Wii_standard;
        % The numbers that sepcify the amount of randomness
        sigma_ee_asym; 
        sigma_ee_rec; 
        sigma_ie; 
    end
    
    %% METHODS THAT DEFINE CONNECTIONS
    methods
        % ------------------------------------------------------------
        function [this]=ConnectionProperties()
            fprintf('Set the connection input properties...\n'); 
        end
        % ------------------------------------------------------------
        function [this,W]=generateConnections(this)
            % Setup shortcut alias
            t=this;
            
            % Create random stream
            r = RandStream.create('mrg32k3a',...
                'NumStreams',1,'Seed','shuffle');
            
            % Get indices of the exitatory and inhibitory cells
            exc         = t.identities == 0;
            inh         = t.identities == 1;
            nNeurons    = numel(t.identities);
            
            % Initialize all of the various weight matrices to be the size
            % of the eventual W matrix
            WEE = zeros(nNeurons,nNeurons);
            WEI = zeros(nNeurons,nNeurons);
            WII = zeros(nNeurons,nNeurons);
            WIE = zeros(nNeurons,nNeurons);
            % Then the eventual recepticle of all these instatiated values
            t.W = zeros(nNeurons,nNeurons);
            
            % Find the diagonal indices for the excitatory and inhibitory
            % cells directly
            exc_diagonal = sub2ind(size(t.W),find(exc),find(exc));
            inh_diagonal = sub2ind(size(t.W),find(inh),find(inh));
            
            %% Compute E to E connections
            % First apply the asymmetric connections to all of the e to e
            % connections
            WEE(exc,exc) = t.Wee_asym + ...
                t.sigma_Wee*(rand(r,1)-0.5); % All the asymmetric different to different receive the same value as a consequence of this
            WEE(exc_diagonal) = t.Wee_recurrent + ...
                t.sigma_ee_recurrent.*rand(r,1,numel(exc_diagonal));
            % Provide all the excitatory to inhibitory connections
            WIE(inh,exc) = t.Wie_standard + ...
                t.sigma_ie(r,size(WIE(inh,exc)));
            WEI(exc,inh) = t.Wei_standard;
            WII(inh,inh) = t.Wii_standard;
            
            
            %% Set total connection matrix
            try
                if gpuDeviceCount > 0
                    t.W = gpuArray(WEE + WIE + WII + WEI);
                else
                    t.W = gpuArray(WEE + WIE + WII + WEI);
                end
            catch
                % if gpuDeviceCount function not found or if gpuArray fails
                % to initialize because of cuda driver problems, then it
                % catches here and simply creates a standard array.
                warning('Caught error while trying to assign output');
                t.W = gpuArray(WEE + WIE + WII + WEI);
            end
            
            % Set the matrix that can be optionally returned
            W=t.W; % for return purposes only
            this=t;
        end
        % ------------------------------------------------------------
    end
    
end