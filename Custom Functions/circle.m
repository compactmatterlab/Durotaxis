function h = circle(x,y,r,varargin)
c = [];
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

for i=1:size(varargin,2)
    if strcmp(varargin{i},'color')
        c = varargin{i+1};
        i=i+1;
    elseif strcmp(varargin{i},'filled')
        fill(xunit, yunit, [1,1,1])
    end
end

h = gobjects;
for i=1:size(xunit,1)
    h(i) = plot(xunit(i,:), yunit(i,:),c);
end
end