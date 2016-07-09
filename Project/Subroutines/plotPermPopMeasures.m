function plotPermPopMeasures(axs, PS, a, b, c, vidFilename)
% Inputs: Number label of two to three conditions to plot against
% eachother. If three are given, a 3D scatter is given instead. But the
% first condition listed is the one by which lines are ordered and plotted.

%% Pre-processing
% 

axisLabels = {'Context','Motion','Color'}; % normally, I don't like hard-coding in things to functions, but time is short

if isempty(axs)
    axs = gca;
end
if nargin < 6
    vidFilename = '3axisVideo.mp4';
end

spec = PS.specifics;
oax = PS.orthogConditionAxes;

% Create matrix of available permutations of trial stimuli
cMat = {spec.cValue}';
cMat = cat(1,cMat{:});

%% Determine mode
if nargin == 4
    plotmode = 1;
elseif nargin >= 5
    plotmode = 2;
end

%% Call 
if plotmode == 1
    TwoDimPlot(a,b);
elseif plotmode == 2
    ThreeDimPlot(a,b,c);
end

%% 2D graph in Time

    function [] = TwoDimPlot(c1,c2)
        [sMat,sortInd1] = sortrows(cMat,c2);
        [sMat,sortInd2] = sortrows(sMat,c1);
        axes(axs);
        
        hold on; legendLab=cell(1,size(sMat,1));
        p = gobjects(size(sMat,1),1);
        for i = 1:size(sMat,1)
            
            proj = spec(sortInd2(sortInd1(i))).proj;
            
            p(i)=plot( proj(:,c1), proj(:,c2) );
            start=plot( proj(1,c1), proj(1,c2), 'color', p(i).Color, 'marker', '*');
            stop=plot( proj(end,c1), proj(end,c2), 'color', p(i).Color,  'marker', '*');
            xlabel(axisLabels{c1});
            ylabel(axisLabels{c2});
            
            % provide arrows to show direction of time, if that function is
            % in path
            if exist('arrowh.h','file') && ~isempty(p(i))
                p_arrow=...
					plotAndCorrectArrows(p(i),proj(:,1),proj(:,2));
            end
            
            legendLab{i} = num2str(sMat(i,[c1 c2]));
            
        end
        
        legend(p,legendLab);
        
    end

%% 3D Scatter in Time

    function [] = ThreeDimPlot(c1,c2,c3)      
        
        [sMat, sortInd]= sort(cMat,c1);
        axes(axs); currf=gcf;
        currf.WindowStyle='normal';
        hold on;
              
        % Initialize animatedline objects
        h = gobjects(1,size(sMat,1));
        for i = 1:size(sMat,1)
            proj = spec(sortInd(i)).proj;
            h(i) = animatedline(proj(1,c1), proj(1,c2), proj(1,c3),...
                'color',rand(1,3));
        end
        xlabel(axisLabels{c1});
        ylabel(axisLabels{c2});
        zlabel(axisLabels{c3});
        view([45 45]);
        forig=getframe(currf);
        x=size(forig.cdata,1);y=size(forig.cdata,2);
        
         v = VideoWriter(vidFilename);
         open(v);
        
        % Move each line
        for t = 1:10:size(spec(1).data(:,1))
            for i = 1:size(sMat,1)
                proj = spec(sortInd(i)).proj;
                addpoints(h(i), proj(t,c1), proj(t,c2), proj(t,c3));
                drawnow;
            end
            camorbit(-pi/120,0);
            f=getframe(currf);
            writeVideo(v,f);
        end
        
        close(v);
        
    end

    %%   plotAndCorrectArrows(p,t,x,y)
	function arrow_handles = plotAndCorrectArrows(p,x,y)
        %
        % Purpose: 
        
        % test if we can throw arrows on the xy data
        if numel(x) < 2 || (numel(x) == 2 && x(1) == x(2) && ...
                y(1) == y(2))
            return;
        end 
        
        traj_color = p.Color; 
        if norm(traj_color,2)~=0 
            traj_color = traj_color/norm(traj_color,2);
        end
        
        arrow_handles=arrowh( x, y, ...
                    traj_color, [80,60], 0:arrowdensity:100);
        
        for ah = 1:numel(arrow_handles)
            obj = handle(arrow_handles(ah));
            obj.FaceAlpha = 0.75;
            obj.EdgeColor = 'k';
            obj.EdgeAlpha = 0.2;
        end
	end
    
end