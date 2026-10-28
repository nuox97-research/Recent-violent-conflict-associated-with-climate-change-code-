function [ varargout ] = heatmapcompute( varargin)
%CLOUDPLOT Does a cloud plot of the data in X and Y.
%
% a modified function based on cloudPlot function by Daniel Armyr 
% (https://www.mathworks.com/matlabcentral/fileexchange/23238-cloudplot)
% by adding two extra arguments for local averaging (true/false, window size) 
% example: heatmapcompute(px,py,pz,datarange,false,[100,100],true,[7,7],1)

[args,pvpairs] = parseparams(varargin);

nArgs = numel(args);
bins = [];                                      % Default value
useLogScale = false;                            % Default value
axisLimitsSet = false;                          % Default value

if ( nArgs >= 1 ); X = args{1}; end             % x factor
if ( nArgs >= 2 ); Y = args{2}; end             % y factor
if ( nArgs >= 3 ); Z = args{3}; end             % z factor
if ( nArgs >= 4 )
    axisLimits = args{4};                       % ranges of x-y axes
    if ( ~isempty(axisLimits) )
        pointSelect = X<=axisLimits(2) & X>=axisLimits(1) & ...
        Y<=axisLimits(4) & Y>=axisLimits(3);
        X = X(pointSelect);
        Y = Y(pointSelect);
        Z = Z(pointSelect);
        axisLimitsSet = true;
    end
end
if ( nArgs >= 5 ); useLogScale = args{5}; end   % use log scale or not
if ( nArgs >= 6 ); isvisible = args{6}; end     % show the image or not
if ( nArgs >= 7 ); bins = args{7}; end          % resolution of output image
if ( nArgs >= 8 ); usesmooth = args{8}; end     % use smooth or not
if ( nArgs >= 9 ); smoothsize = args{9}; end    % smoothing window size 


assert ( nArgs < 10, 'No more than 9 arguments plus name-value pairs allowed' );
assert ( nArgs >= 2, 'At least 2 arguments required' );

% Check the data size
assert ( numel(X) == numel(Y), ...
    'The number of elements in X and Y must be the same.' );
assert ( numel(Y) == numel(Z), ...
    'The number of elements in Y and Z must be the same.' );
assert ( numel(X) == numel(Z), ...
    'The number of elements in X and Z must be the same.' );
  
%Remove any nans or Infs in the data as they have no meaning in this
%context.
pointSelect = ~(isinf(X) | isnan(X) | isinf(Y) | isnan(Y) | isinf(Z) | isnan(Z) );
X = X(pointSelect);
Y = Y(pointSelect);
Z = Z(pointSelect);

% Plot to get appropriate limits
h = [];
if ( axisLimitsSet )    
    g = gca;
    set ( g, 'Xlim', [axisLimits(1) axisLimits(2)] );
    set ( g, 'Ylim', [axisLimits(3) axisLimits(4)] );
    set ( g, 'units', 'pixels' );
else
    h = plot ( X(:), Y(:), '.' );        
    g = get( h, 'Parent' );        
end
xLim = get(g, 'Xlim' );
yLim = get(g, 'Ylim' );

%Get the bin size.
unitType = get(g,'Unit');
set(g,'Unit','Pixels')
axesPos = get(g,'Position');
nHorizontalBins = axesPos(3);
nVerticalBins = axesPos(4);
set(g,'Unit', unitType );

% Clear the data, as we actually don't want to see it.
if ( ~isempty(h) )
    set ( h, 'XData', [] );
    set ( h, 'YData', [] );
end

% Allocate an area to draw on
if ( isempty(bins) )
    bins = ceil([nHorizontalBins nVerticalBins ]);
else
    assert ( isnumeric(bins) && isreal(bins), ...
        'Bin counts must be real numbers' );
    assert ( numel(bins) == 2, ...
        'There must be exactly 2 bin values, one for x and one for y' );
    if( bins(1) > nHorizontalBins && bins(2) > nVerticalBins )
        warning ( 'CLOUDPLOT:TOO_MANY_BINS', ...
            'Having bins smaller than one pixel on the screen is not advised' );
    end        
        
end
binSize(2) = diff(yLim)./(bins(2));
binSize(1) = diff(xLim)./(bins(1));

numpoint = zeros(bins);
denpoint = zeros(bins);
% Draw in the canvas
xBinIndex = floor((X - xLim(1))/binSize(1))+1;
yBinIndex = floor((Y - yLim(1))/binSize(2))+1;

% Added security: Make sure indexes are never outside canvas. May not be
% possible.
pointsSelect = xBinIndex > 0 & xBinIndex <= bins(1) & ...
    yBinIndex > 0 & yBinIndex <= bins(2);
xBinIndex = xBinIndex(pointsSelect);
yBinIndex = yBinIndex(pointsSelect);

for i = 1:numel(xBinIndex)
    numpoint(xBinIndex(i),yBinIndex(i)) = ...
        numpoint(xBinIndex(i),yBinIndex(i)) + 1;
    denpoint(xBinIndex(i),yBinIndex(i)) = ...
        denpoint(xBinIndex(i),yBinIndex(i)) + Z(i);
end

denpoint = denpoint./numpoint;
denpoint( isinf(denpoint) ) = 0;  %%%%%%%%%%%%%%%%% or nan

% Show the canvas and adjust the grids.
if ( useLogScale )
    if (usesmooth)
       Hfilter = fspecial('average',smoothsize);
       denpoint = imfilter(denpoint,Hfilter);        
    end
    if (isvisible) 
        h = imagesc(xLim, yLim,log10(denpoint)');
        set(h,'alphadata',~isnan(denpoint'));
    end
else
    if (usesmooth)
       Hfilter = fspecial('average',smoothsize);
       denpoint = imfilter(denpoint,Hfilter);        
    end
    if (isvisible)
        h = imagesc(xLim, yLim, denpoint');    
        set(h,'alphadata',~isnan(denpoint'));
    end
end

axis ( 'xy' );
axis ( 'tight' );
set ( g, 'units', 'normalized' ); % Enable resizing 

% Optionally return a handle.
if ( nargout >= 1)
    varargout{1} = h;
end
if ( nargout >= 2)
    varargout{2} = denpoint;
end
