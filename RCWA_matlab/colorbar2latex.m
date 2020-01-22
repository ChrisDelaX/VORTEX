%%% by Christian Delacroix (2014)



function colorbar2latex(hC)


if strcmp(get(hC,'XDir'),'reverse')
set(hC,'XDir','normal')
xReverseFlag = 1; 
else
xReverseFlag = 0;
end
if strcmp(get(hC,'YDir'),'reverse')
set(hC,'YDir','normal')
yReverseFlag = 1; 
else
yReverseFlag = 0;
end
xLabPosCurr = get(get(hC,'XLabel'),'Position');
yLabPosCurr = get(get(hC,'YLabel'),'Position');


xLabelDy = 0;
yLabelDx = 0;
hAxis = hC;


% Get properties for the specified axis
xLimit = get(hAxis,'XLim');
xTick = get(hAxis,'XTick');
xTickLabel = get(hAxis,'XTickLabel');
xTickLabelMode = get(hAxis,'XTickLabelMode');
xScale = get(hAxis,'XScale');
xAxisLocation = get(hAxis,'XAxisLocation');

yLimit = get(hAxis,'YLim');
yTick = get(hAxis,'YTick');
yTickLabel = get(hAxis,'YTickLabel');
yTickLabelMode = get(hAxis,'YTickLabelMode');
yScale = get(hAxis,'YScale');
yAxisLocation = get(hAxis,'YAxisLocation');

if ~exist('fontSize','var')
    fontSize = get(hAxis, 'FontSize');
end

% Get properties for the current axis
xLimitCurr = get(gca,'XLim');
xScaleCurr = get(gca,'XScale');

yLimitCurr = get(gca,'YLim');
yScaleCurr = get(gca,'YScale');

% Clear the current labels
set(hAxis,'XTickLabel',[]);
set(hAxis,'YTickLabel',[]);

% Get position of the specified axis
posAxis = get(hAxis,'Position');

% Get position of the current axis
posCurr = get(gca, 'Position');

% Convert x point between figure and axis data coordinates on linear scale
xFig2DatLinAxis = @(x)(xLimit(1) + diff(xLimit)*(x-posAxis(1))/posAxis(3));
xDat2FigLinAxis = @(x)(posAxis(1) + (x - xLimit(1))*posAxis(3)/diff(xLimit));

xFig2DatLinCurr = @(x)(xLimitCurr(1) + diff(xLimitCurr)*(x-posCurr(1))/posCurr(3));
xDat2FigLinCurr = @(x)(posCurr(1) + (x - xLimitCurr(1))*posCurr(3)/diff(xLimitCurr));

% Convert y point between figure and axis data coordinates on linear scale
yFig2DatLinAxis = @(y)(yLimit(1) + diff(yLimit)*(y-posAxis(2))/posAxis(4));
yDat2FigLinAxis = @(y)(posAxis(2) + (y - yLimit(1))*posAxis(4)/diff(yLimit));

yFig2DatLinCurr = @(y)(yLimitCurr(1) + diff(yLimitCurr)*(y-posCurr(2))/posCurr(4));
yDat2FigLinCurr = @(y)(posCurr(2) + (y - yLimitCurr(1))*posCurr(4)/diff(yLimitCurr));

% Convert x point between figure and axis data coordinates on logarithmic scale
xFig2DatLogAxis = @(x)(exp(log(xLimit(1)) + log(xLimit(2)/xLimit(1))*(x-posAxis(1))/posAxis(3)));
xDat2FigLogAxis = @(x)(posAxis(1) + posAxis(3)*log(x/xLimit(1))/log(xLimit(2)/xLimit(1)));

xFig2DatLogCurr = @(x)(exp(log(xLimitCurr(1)) + log(xLimitCurr(2)/xLimitCurr(1))*(x-posCurr(1))/posCurr(3)));
xDat2FigLogCurr = @(x)(posCurr(1) + posCurr(3)*log(x/xLimitCurr(1))/log(xLimitCurr(2)/xLimitCurr(1)));

% Convert y point between figure and axis data coordinates on logarithmic scale
yFig2DatLogAxis = @(y)(exp(log(yLimit(1)) + log(yLimit(2)/yLimit(1))*(y-posAxis(2))/posAxis(4)));
yDat2FigLogAxis = @(y)(posAxis(2) + posAxis(4)*log(y/yLimit(1))/log(yLimit(2)/yLimit(1)));

yFig2DatLogCurr = @(y)(exp(log(yLimitCurr(1)) + log(yLimitCurr(2)/yLimitCurr(1))*(y-posCurr(2))/posCurr(4)));
yDat2FigLogCurr = @(y)(posCurr(2) + posCurr(4)*log(y/yLimitCurr(1))/log(yLimitCurr(2)/yLimitCurr(1)));

% Convert x point between figure and axis [0,1] coordinates
xFig2Ax = @(x)((x - posAxis(1))/posAxis(3));
xAx2Fig = @(x)(posAxis(1) + x*posAxis(3));

% Convert y point between figure and axis [0,1] coordinates
yFig2Ax = @(y)((y - posAxis(2))/posAxis(4));
yAx2Fig = @(y)(posAxis(2) + x*posAxis(4));

switch xScale
    case 'linear'
        xFig2DatAxis = xFig2DatLinAxis;
        xDat2FigAxis = xDat2FigLinAxis;
    case 'log'
        xFig2DatAxis = xFig2DatLogAxis;
        xDat2FigAxis = xDat2FigLogAxis;
end

switch yScale
    case 'linear'
        yFig2DatAxis = yFig2DatLinAxis;
        yDat2FigAxis = yDat2FigLinAxis;
    case 'log'
        yFig2DatAxis = yFig2DatLogAxis;
        yDat2FigAxis = yDat2FigLogAxis;
end

switch xScaleCurr
    case 'linear'
        xFig2DatCurr = xFig2DatLinCurr;
        xDat2FigCurr = xDat2FigLinCurr;
    case 'log'
        xFig2DatCurr = xFig2DatLogCurr;
        xDat2FigCurr = xDat2FigLogCurr;
end

switch yScaleCurr
    case 'linear'
        yFig2DatCurr = yFig2DatLinCurr;
        yDat2FigCurr = yDat2FigLinCurr;
    case 'log'
        yFig2DatCurr = yFig2DatLogCurr;
        yDat2FigCurr = yDat2FigLogCurr;
end


if ~isempty(yTickLabel)
    % Set the Y axis
    
    yTickIndex = find((yTick >= yLimit(1)) & (yTick <= yLimit(2)));
    yTickVisible = yTick((yTick >= yLimit(1)) & (yTick <= yLimit(2)));
    yLabel = cell(length(yTickVisible),1);
    
    switch yTickLabelMode
        case 'auto'
            assert(length(yTickVisible) <= size(yTickLabel,1));
            switch yScale
                case 'linear'
                    % Determine where there should be an exponent
                    yExp = max(abs(yLimit));
                    if yExp > 0
                        yExpLog = ceil(log10(yExp))-1;
                        if(yExpLog > 0)
                            yExpLog = yExpLog - 1;
                        end
                        if (yExpLog > -3) && (yExpLog <= 3)
                            yExpLog = 0;
                        end
                    else
                        yExpLog = 0;
                    end
                    yExp = 10^yExpLog;
                    
                    for i=1:length(yTickVisible)
                        value = yTickVisible(i)/yExp;
                        if abs(value) <= eps
                            value = 0;
                        end
                        yLabel{i} = ['$' num2str(value) '$'];
                    end
                    
                    if (abs(yExpLog) > eps) && (abs(yExpLog) < 1/eps)
                        switch yAxisLocation
                            case 'left'
                                if xReverseFlag==0
                                    hText = text(...
                                        xFig2DatCurr(xDat2FigAxis(xLimit(1))),...
                                        yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03),...
                                        ['$\times 10^{' num2str(yExpLog) '}$'],...
                                        'Interpreter','latex',...
                                        'FontSize', fontSize);
                                elseif xReverseFlag==1 % XDir reverse
                                    hText = text(...
                                        xFig2DatCurr(xDat2FigAxis(xLimit(2))),...
                                        yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03),...
                                        ['$\times 10^{' num2str(yExpLog) '}$'],...
                                        'Interpreter','latex',...
                                        'FontSize', fontSize);
                                end
                            case 'right'
                                if xReverseFlag==0
                                    hText = text(...
                                        xFig2DatCurr(xDat2FigAxis(xLimit(2))),...
                                        yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03)+10,...
                                        ['$\times 10^{' num2str(yExpLog) '}$'],...
                                        'HorizontalAlignment','Right',...
                                        'Interpreter','latex',...
                                        'FontSize', fontSize);
                                elseif xReverseFlag==1 % XDir reverse
                                    hText = text(...
                                        xFig2DatCurr(xDat2FigAxis(xLimit(1))),...
                                        yFig2DatCurr(yDat2FigAxis(yLimit(2)) + 0.03),...
                                        ['$\times 10^{' num2str(yExpLog) '}$'],...
                                        'HorizontalAlignment','Right',...
                                        'Interpreter','latex',...
                                        'FontSize', fontSize);
                                end
                        end
                        set(hText,'HitTest','off');
                    end
                case 'log'
                    for i=1:length(yTickVisible)
                        sgn = sign(yTickVisible(i));
                        yExp = log10(yTickVisible(i));
                        if abs(yExp) <= eps
                            yExp = 0;
                        end
                        yLabel{i} = ['$' num2str(10*sgn) '^{' num2str(yExp) '}$'];
                    end
            end
        case 'manual'
            if iscell(yTickLabel)
                for i=1:length(yTickVisible)
                    yLabel{i} = yTickLabel{1+mod(yTickIndex(i)-1,length(yTickLabel))};
                end
            else
                for i=1:length(yTickVisible)
                    yLabel{i} = yTickLabel(1+mod(yTickIndex(i)-1,size(yTickLabel,1)),:);
                end
            end
    end
    
    xMin = bitmax;
    
    switch yAxisLocation
        case 'left'
            for i = 1:length(yTickVisible)
                if xReverseFlag==0
                    hText = text(...
                        xFig2DatCurr(xDat2FigAxis(xLimit(1))-0.01),...
                        yFig2DatCurr(yDat2FigAxis(yTickVisible(i))),...
                        strtrim(yLabel{i}),...
                        'HorizontalAlignment','Right',...
                        'Interpreter','latex',...
                        'FontSize', fontSize);
                elseif xReverseFlag==1 % XDir reverse
                    hText = text(...
                        xFig2DatCurr(xDat2FigAxis(xLimit(2))+0.01),...
                        yFig2DatCurr(yDat2FigAxis(yTickVisible(i))),...
                        strtrim(yLabel{i}),...
                        'HorizontalAlignment','Right',...
                        'Interpreter','latex',...
                        'FontSize', fontSize);
                end
                set(hText,'HitTest','off');
                set(hText,'Units','normalized');
                xExt = get(hText,'Extent');
                set(hText,'Units','data');
                xMin = min(xMin, xExt(1));
            end
        case 'right'
            for i = 1:length(yTickVisible)
                if xReverseFlag==0
                    hText = text(...
                        xFig2DatCurr(xDat2FigAxis(xLimit(2))+0.005),...
                        yFig2DatCurr(yDat2FigAxis(yTickVisible(i))+0.005),...
                        strtrim(yLabel{i}),...
                        'HorizontalAlignment','Left',...
                        'Interpreter','latex',...
                        'FontSize', fontSize);
                elseif xReverseFlag==1 % XDir reverse
                    hText = text(...
                        xFig2DatCurr(xDat2FigAxis(xLimit(1))-0.01),...
                        yFig2DatCurr(yDat2FigAxis(yTickVisible(i))),...
                        strtrim(yLabel{i}),...
                        'HorizontalAlignment','Left',...
                        'Interpreter','latex',...
                        'FontSize', fontSize);
                end
                set(hText,'HitTest','off');
                set(hText,'Units','normalized');
                xExt = get(hText,'Extent');
                set(hText,'Units','data');
                xMin = min(xMin, xExt(1));
            end
    end
    
    yLabel = get(hAxis,'YLabel');
    yLabelPos = get(yLabel,'Position');
    set(yLabel,'Position',[xFig2DatCurr(xAx2Fig(xMin) - 0.002 + yLabelDx)+10 yLabelPos(2) yLabelPos(3)])
    
    ylim(hAxis, yLimit);
end


% reverse axes
if xReverseFlag
    set(hC,'XDir','reverse')
    set(get(hC,'ylabel'),'position',[2*xLimitCurr(1)+diff(xLimitCurr)-yLabPosCurr(1) 2*yLimitCurr(1)+diff(yLimitCurr)-yLabPosCurr(2) yLabPosCurr(3)])
end
if yReverseFlag
    set(hC,'YDir','reverse')
    set(get(hC,'xlabel'),'position',[2*xLimitCurr(1)+diff(xLimitCurr)-xLabPosCurr(1) 2*yLimitCurr(1)+diff(yLimitCurr)-xLabPosCurr(2) xLabPosCurr(3)])
end


end










