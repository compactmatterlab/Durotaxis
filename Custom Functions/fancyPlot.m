

function pl = fancyPlot(varargin) %(xdata, ydata, color | title | xlabel | ylabel | legend | lineStyle | marker | errorBars | lineWidth | fontSize)

plotColors = 0;
setLegend = 0;
setLineSpec = 0;
setMarker = 0;
getSpline = 0;
setErrorBars = 0;
lnWidth = 3;
fntSize = 24;

%figure
%set(gcf, 'Position', [20, 50, 900, 700])
hold on

if size(varargin,2)>2
    for i=3:size(varargin,2)
        str = varargin{i}{1};
        if strcmp(str,'color')
            plotColors = 1;
            for j=2:size(varargin{i},2)
                col(j-1,:) = varargin{i}{j}; 
            end            
        elseif strcmp(str,'lineStyle')
            setLineSpec =1;
            for j=2:size(varargin{i},2)
                lSpec(j-1) = {varargin{i}{j}}; 
            end 
        elseif strcmp(str,'marker')
            setMarker = 1;
            for j=2:size(varargin{i},2)
                mark(j-1) = {varargin{i}{j}}; 
            end 
        elseif strcmp(str,'markerSize')
            setMSize = 1;
            for j=2:size(varargin{i},2)
                mSize(j-1) = varargin{i}{j}; 
            end 
        elseif strcmp(str,'title')
            ttl = varargin{i}{2};
            title(ttl)
        elseif strcmp(str,'xlabel')
            xlab = varargin{i}{2};
            xlabel(xlab)
        elseif strcmp(str,'xlim')
            xl = varargin{i}{2};
            xlim(xl)
        elseif strcmp(str,'ylim')
            yl = varargin{i}{2};
            ylim(yl)
        elseif strcmp(str,'ylabel')
            ylab = varargin{i}{2};
            ylabel(ylab)
        elseif strcmp(str,'legend')
            setLegend = 1;
            for j=2:size(varargin{i},2) 
                legendList{j-1} = varargin{i}{j};
            end
        elseif strcmp(str,'spline')
            getSpline = 1;
        elseif strcmp(str,'errorBars')
            setErrorBars = 1;
            for j=2:size(varargin{i},2)
                errbars(j-1) = {varargin{i}{j}}; 
            end  
        elseif strcmp(str,'lineWidth')
            lnWidth = varargin{i}{2};
        elseif strcmp(str,'fontSize')
            fntSize = varargin{i}{2};
        end
    end
end
for i=1:size(varargin{1},2)
    clear xdata ydata
    xdata = varargin{1}{i};
    if ~iscolumn(xdata)
        xdata = xdata';
    end
    ydata = varargin{2}{i};
    if ~iscolumn(ydata)
        ydata = ydata';
    end

    if getSpline
        x2 = min(xdata):(max(xdata)-min(xdata))/100:max(xdata);
        [f, ~] = fit( xdata, ydata, 'pchipinterp', 'Normalize', 'on' );
        ysp = f(x2);
        xsp = x2;
        p = plot(xsp,ysp,'LineWidth',lnWidth);
        try
            for j=1:length(xdata)
                ind(j) = find(abs(xsp-xdata(j))<1e-12);
            end
            p.MarkerIndices = ind;
            
        catch
            p.MarkerIndices = round(linspace(1,length(xsp),length(xdata)));
        end
         %disp(p.MarkerIndices)
    else
        p = plot(xdata,ydata,'LineWidth',lnWidth);
    end
    
    pl(i) = p;
    xl = xlim;
    yl = ylim;
    


    if plotColors
        p.Color = col(i,:);
    end
    if setLineSpec
        p.LineStyle = lSpec{i};
    end
    if setMarker
        p.Marker = mark{i};
        if plotColors
            p.MarkerFaceColor = col(i,:);
        else
            p.MarkerFaceColor = p.Color;
        end
        if setMSize
            p.MarkerSize = mSize(i);
        end
    end
    
    if setErrorBars
        e = errorbar(xdata,ydata,errbars{i});
        e.LineStyle = 'none';
        e.Color = p.Color;
        e.MarkerFaceColor = p.Color;
        e.MarkerEdgeColor = p.Color;
        e.MarkerSize = p.MarkerSize;
        e.CapSize = lnWidth*4;
        e.LineWidth = lnWidth;
    end
    
end
if setLegend
    l=legend(pl(1:size(varargin{1},2)),legendList(1:size(varargin{1},2)),'Location','best');
    l.EdgeColor = [1,1,1];
    l.FontSize = fntSize;
    set(l,'EdgeColor','none');
    set(l,'color','none');
end
box on
set(gca,'fontsize',fntSize,'FontName', 'Calibri')
set(gca,'color','none')


end