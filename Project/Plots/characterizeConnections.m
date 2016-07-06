function characterizeConnections(W,neurIdentities)
    % Generate plot of connections
            out=figure(100); clf;
            subplot(1,2,1);
            axs = 1:size(W,1)
            laxs=[axs(1)-0.5 axs axs(end)+0.5];
            imagesc(axs,axs,W);
            colorbar;
            title('Connectivity Matrix');
            xlabel('Neuron #');
            ylabel('Neuron #');
            % Demarcate EI border, if neurIdentites were grouped by E and
            % I.
            if nargin == 2
                loc = diff(neurIdentities);
                if sum(loc) == 1
                    loc = find(loc)+0.5;
                    hold on;
                    l=line(repmat(loc,size(laxs)),laxs);
                    set(l,'color','r','linewidth',2,'linestyle',':');
                    l=line(laxs,repmat(loc,size(laxs)));
                    set(l,'color','r','linewidth',2,'linestyle',':');
                end
            end
                
            subplot(1,2,2);
            laxs = 1:size(W,1)
            logmeasure=log10(abs(W));
            logmeasure(logmeasure==-inf) = 0;
            logmeasure = logmeasure.*sign(W);
            imagesc(axs,axs,logmeasure);
            colorbar;
            title('''SignedLog'' Connectivity Measure');
            xlabel('Neuron #');
            ylabel('Neuron #');
            % Demarcate EI border, if neurIdentites were grouped by E and
            % I.
            if nargin == 2
                loc = diff(neurIdentities) ~=0;
                if sum(loc) == 1
                    loc = find(loc)+0.5;
                    hold on;
                    l=line(repmat(loc,size(laxs)),laxs);
                    set(l,'color','r','linewidth',2,'linestyle',':');
                    l=line(laxs,repmat(loc,size(laxs)));
                    set(l,'color','r','linewidth',2,'linestyle',':');
                end
            end
end