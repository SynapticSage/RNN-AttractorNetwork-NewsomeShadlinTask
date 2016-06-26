function fillHandle = plot_digital(mandArg,varargin)
% A function that plots a digital signal as a shaded region onto an axis.
% Digital is defined as zero versus non-zero signal. This is something I do
% fairly often and figured that it was time to make an official function
% for this type of plot.
%
% TODO :    Make sure that this can handle a bunch of signals provided in the
%           columns all at once.

%% Decide how to handle arguments

% Set defaults
axs=gca;
taxis=[];
colorspec = 'k';
transp = 0.25;

% Parse the arguments
if nargin > 0
    if nargin == 1
        
        signal = mandArg;
        
    elseif narargin == 2
        
        taxis = mandArg;
        signal = varargin{1};
        
    else % MORE THAN 2 ARGUMENTS (I.E. there're optional args)
        % CONDITION 1: Time and Signal Passed
        if ~ischar(varargin{1})
            taxis = mandArg;
            signal=varargin{1};
            varargin=varargin(2:end);
        % CONDITION 2: Only Signal Passed
        else
            signal = mandArg;
        end
        
        % Parse optional arguments
        for v = 1:2:numel(varargin)
            switch varargin{v}
                case 'axis',        axs = varargin{v+1};
                case 'transp',      transp = varargin{v+1};
                case 'colorspec',   colorspec = varargin{v+1};
                otherwise, error('Unrecognized input!');
            end
        end
    end
else
    error('Inputs wrong!');
end
    
% If time axis wasn't given, then set it the index
if isempty(taxis), taxis = 1:numel(signal); end;

% detect if signal is a row, and if it is, orient it to a column to be
% consistent with the code below
if isrow(signal)
    signal = signal';
end
if isrow(taxis)
    taxis = taxis';
end

%% Generate plot

% Store the Ylimits of the axis...
YLim = get(axs,'YLim');

% Now we need to find all times that the signal is on
onTimes = signal ~= 0;

% Now generate high and low segments
upper = onTimes * YLim(2);
lower = onTimes * YLim(1);

% Reverse the direction of the signal for the lower, we have to create a
% closed region in order to make the fill option work
lower = lower(end:-1:1);

% Now create the polygon we have hto fill
Y = [upper; lower];
X = [taxis; taxis(end:-1:1)];

% Fill the area
axis(axs);
fillHandle = fill(X,Y,colorspec);
fillHandle.EdgeAlpha = transp;
fillHandle.FaceAlpha = transp;


end