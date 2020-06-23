clear global
clearvars ~except
close all 
global m b dx rto fldrOpen prefix
load('GNetData2')
load('GNetData3')




%----------------------------------------------------------------------------------
%Values you need to set
%----------------------------------------------------------------------------------
%Set this each time for the video you want to analyze
fileNum = 1;
autoUpdate = 0;
 
%Change these for the cell type and date of imaging
prefix = 'MDA-MB231_Bleb_Jan92020';
oldName = prefix;

%Change this file path for your computer
%fldrSave = sprintf('C:/Users/Ben/Google Drive/Research/%s Results/%s', prefix(1:find(prefix=='_')-1),prefix(find(prefix=='_')+1:end));
fldrSave = sprintf('Z:/Ben Yeoman/Timelapse Images/Bleb/%s',prefix(find(prefix=='_')+1:end));
fldrOpen = sprintf('Z:/Ben Yeoman/Timelapse Images/Bleb/%s',prefix);

%Automatically upates filenumber for next video
str = sprintf('%s/S%dData',fldrSave,fileNum);
if autoUpdate
    while isfile([str, '.mat'])
        fileNum=fileNum+1;
        str = sprintf('%s/S%dData',fldrSave,fileNum);
    end
end
disp(['Using ', str])

%Organize files if needed
if numel(dir([fldrOpen,sprintf('/S%d',fileNum)]))<3
    Organize(oldName,240)
end



%----------------------------------------------------------------------------------

%----------------------------------------------------------------------------------
% Commands
% Left click - Add new cell position
% Right click - Move existing cell position
% Space - Goto next image
% Backspace goto previous 
% S - Inverts cell selction window image (Optimizes auto-select for if cell 
% are either bright-centered or dark-centered
% D - Lock closest cell position
% F - Unlock closest cell position
% H - Hide cell tracks
% J - Show cell tracks
% Del - Delete closest cell
%----------------------------------------------------------------------------------


width = 60;
rad = 30;
r = 15;
xres = 512;
yres = 512;
if xres == 512
    len = 88;
    rto = 2.55;
else
    len = 190;
    rto = 2.6;
end


play = 1;



fileIn = sprintf('%s/S%d/%s_s%d_t',fldrOpen,fileNum,prefix,fileNum);
fileOut = sprintf('%s/S%d/%s_v%d_t',fldrOpen,fileNum,prefix,fileNum);


if exist(fldrSave,'dir') ~= 7 
    mkdir(fldrSave)
end

O = zeros(yres,xres,97);
% K = zeros(yres,xres,3,97);
% Q = zeros(yres,xres,3,97);
J = zeros(yres,xres);
P = zeros(yres,xres);
M = zeros(yres,xres);
% R = zeros(yres,xres);
n = 0;


%Define plot colors
for i=1:5
    col((i-1)*7+1:(i-1)*7+7,:) = get(gca,'ColorOrder');
end

for i=1:97
    str = sprintf('%s%d.tif',fileIn,i);
    O(:,:,i) = rescale(imread(str));
end
if play
    ip=implay(O,8);
    if play==1
        ip.Parent.Position=[1 41 1536 724];
    else
        ip.Parent.Position=[1921 -175 1920 964.8];
    end
end
set(gcf, 'Position', [1 41 1536 748.8000]);
im = imshow(O(:,:,1),'InitialMagnification','fit');
hold on

%Gets boudaries from user input
i = 0;
u = zeros;
v = zeros;
t = text(3,506,'Draw Boundaries','Color','white');
s = gobjects;
l1 = gobjects;
l2 = gobjects;
while 1
    i=i+1;
    [u(i,1),v(i,1),btn] = ginput(1);
    if isempty(btn)
        btn = -1;
    end
    s(i) = scatter(u(i,1),v(i,1),10,'w','filled');
    if i > 1 && btn ~= 127 && btn ~= 32 
        xf = [u(end) u(end-1)];
        yf = [v(end) v(end-1)];
        c = [[1; 1]  xf(:)]\yf(:);                        
        m = c(2);
        b = c(1);
        o = atand(1/abs(m));
        dx = len/sind(o);
        b = b-dx;
        lx1(1:5,1) = 0;
        lx1(1:5,2) = xres;
        lx2 = lx1;
        ly1(1:5,1) = b+linspace(-1,3,5)*dx*rto;
        ly1(1:5,2) = m*xres+b+linspace(-1,3,5)*dx*rto;
        ly2 = ly1+dx;
        for j=1:5
            l1(j) = line([lx1(j,1),lx1(j,2)],[ly1(j,1),ly1(j,2)],'Color','w');
            l2(j) = line([lx2(j,1),lx2(j,2)],[ly2(j,1),ly2(j,2)],'Color','w');
        end
    end
    if btn == 127
        delete([s,l1,l2])
        i = 0;
        u = 0;
        v = 0; 
    end
    if btn == 32
        delete([s,t])
        t = text(3,506,'Pick Initial Positions','Color','white');
        break;
    end
end

%Save Boundary Edges
str = sprintf('%s/L%dData',fldrSave,fileNum);
save(str,'lx1','lx2','ly1','ly2')


for i=1:xres
    yt(1:5) = m*i+b+linspace(-1,3,5)*dx*rto;
    yt(6:10) = yt(1:5)+dx;
    for j=1:yres
        if (j > yt(1) && j < yt(6)) || (j > yt(2) && j < yt(7)) || (j > yt(3) && j < yt(8)) || (j > yt(4) && j < yt(9)) || (j > yt(5) && j < yt(10))
            M(j,i) = 1;
        else
            M(j,i) = 0;
        end
    end
end


p = gobjects;
l3 = gobjects;
cellPos=zeros(1,4,1);
div = zeros;
st = 0;
bold = 0;

%Start from previously saved position if available 
try
    load(sprintf('%s/S%dData',fldrSave,fileNum)) 
    div = zeros(size(cellPos,3),1);
    i=size(cellPos,1);
    cellPos(:,:,~cellPos(end,1,:)&~cellPos(end,2,:))=[];
    x(:,1)=cellPos(end-1,1,:);
    y(:,1)=cellPos(end-1,2,:);
    div(:,1)=cellPos(1,4,:);
    for j=1:size(cellPos,3)
        if cellPos(end,4,j) == 0
            cellPos(:,:,j) = [];
            div(j)=[];
            x(j) = [];
            y(j) = [];
        end
    end
catch
    i=1;
    x=zeros;
    y=zeros;
end
while i <= 97

    delete(im);
    im = imshow(O(:,:,i),'InitialMagnification','fit');
    uistack(im,'bottom')
    
    if i > 1 && btn~=8
        clear midx dmin
        dmin = zeros;
        
        P = zeros(yres,xres);
        J = zeros(yres,xres);
        for n=1:numel(x)
            
            for j=1:numel(x)
                dmin(j)=pdist([y(n),x(n);y(j) x(j)]);
            end
            dmin(dmin==0) = [];
            
            %Image Region for Neural Network 
            if x(n) < 1; x(n) = 1;elseif x(n) > xres; x(n) = xres; end 
            if y(n) < 1; y(n) = 1; elseif y(n) > yres; y(n) = yres; end
            if x(n)-rad < 1; x1=1; x2=width; xp = 1-x(n);
            elseif x(n)+rad > xres; x1=xres-width+1; x2=xres; xp = x1+(xres-x(n));
            else; x1=x(n)-rad; x2=x(n)+rad-1; xp = x1;
            end 
            if y(n)-rad < 1; y1=1; y2=width; yp = 1-y(n);
            elseif y(n)+rad > yres; y1=yres-width+1; y2=yres; yp = y1+(yres-y(n));
            else; y1=y(n)-rad; y2=y(n)+rad-1; yp=y1;
            end 
            bd = [y1;y2;x1;x2];
            
            if i==97
                clear img
                img(:,:,1) = O(bd(1):bd(2),bd(3):bd(4),i);
                img(:,:,2) = O(bd(1):bd(2),bd(3):bd(4),i-1);
                pos = predict(net2, img);
                pos = round(pos);
                x(n) = pos(1)+xp;
                y(n) = pos(2)+yp; 
            else
                img(:,:,1) = O(bd(1):bd(2),bd(3):bd(4),i-1);
                img(:,:,2) = O(bd(1):bd(2),bd(3):bd(4),i);
                img(:,:,3) = O(bd(1):bd(2),bd(3):bd(4),i+1);
                pos = predict(net3, img);
                pos = round(pos); 
                x(n) = pos(1)+xp;
                y(n) = pos(2)+yp;  
            end
             
%             %Prepare Local Image Region
%             if x(n) < 1; x(n) = 1; elseif x(n) > xres; x(n) = xres; end 
%             if y(n) < 1; y(n) = 1; elseif y(n) > yres; y(n) = yres; end
%             xprev = x(n);
%             yprev = y(n);
%             if x(n)-rad < 1; x1=1; x2=x(n)+rad; 
%             elseif x(n)+rad > xres; x1=x(n)-rad; x2=xres; 
%             else; x1=x(n)-rad; x2=x(n)+rad; 
%             end 
%             if y(n)-rad < 1; y1=1; y2=y(n)+rad; 
%             elseif y(n)+rad > yres; y1=y(n)-rad; y2=yres; 
%             else; y1=y(n)-rad; y2=y(n)+rad; 
%             end 
%             bd = [y1;y2;x1;x2];
%             I = zeros(bd(2)-bd(1)+1,bd(4)-bd(3)+1);
%             R = zeros(yres,xres);
%             R(y(n),x(n)) = 1;
%             R = imdilate(R,strel('disk',rad,8));      
%             R = R(bd(1):bd(2),bd(3):bd(4)) == 1;
%             edge = bwmorph(R,'remove');
%             if st
%                 It = 1-O(bd(1):bd(2),bd(3):bd(4),i);
%             else
%                 It = O(bd(1):bd(2),bd(3):bd(4),i);
%             end
%             B = imopen(It,strel('disk',15));
%             It = It-B;
%             It(~R) = mean(It(edge));
%             idx = abs(It - mean(It(edge))) < 5*std(It(edge));
%             It(idx) = mean(It(edge));
% 
%                 
%             %Process Local Image Region
%             if min(dmin) > rad
%                 I(:,:,1) = It;
%                 I(:,:,2) = imbinarize(I(:,:,1));
%                 CC = bwconncomp(I(:,:,2));
%                 S = regionprops(CC, 'Area');
%                 if CC.NumObjects > 2
%                     area = round(mean([S.Area]));
%                 else
%                     area = 0;
%                 end
%                 I(:,:,3) = bwareaopen(I(:,:,2), area); 
%                 I(:,:,4) = bwmorph(I(:,:,3),'spur',Inf);
%                 I(:,:,5) = imclose(I(:,:,4),strel('disk',10,8));
%                 I(:,:,6) = imfill(I(:,:,5),'holes');
%                 J(bd(1):bd(2),bd(3):bd(4)) = I(:,:,end);
% 
%                 %Calculate Cell Position
%                 c = regionprops(I(:,:,end),'centroid');
%                 cen = cat(1, c.Centroid);
%                 try
%                     x(n) = round(cen(1))+xprev-rad;
%                     if x(n) < 1; x(n) = 1; elseif x(n) > xres; x(n) = xres; end
%                     y(n) = round(cen(2))+yprev-rad; 
%                     if y(n) < 1; y(n) = 1; elseif y(n) > yres; y(n) = yres; end
%                 catch
%                     x(n) = xprev;
%                     y(n) = yprev;
%                 end     
%             else
%                 clear u v sim
%                 It = O(bd(1):bd(2),bd(3):bd(4),i);
%                 Ip = O(bd(1):bd(2),bd(3):bd(4),i-1);
%                 Ip = Ip(floor(size(Ip,1)/2)-r+1:floor(size(Ip,1)/2)+r+1,floor(size(Ip,2)/2)-r+1:floor(size(Ip,2)/2)+r+1);
%                 R2 = imerode(R,strel('disk',r,8));
%                 [u,v] = ind2sub(size(R2), find(R2));
%                 I = zeros(2*r+1,2*r+1);
%                 for j=1:numel(u)
%                     if v(j)-r < 1; v(j)=r+1; end
%                     if v(j)+r > size(It,1); v(j)=size(It,1)-r; end
%                     if u(j)-r < 1; u(j)=r+1; end
%                     if u(j)+r > size(It,2); u(j)=size(It,2)-r; end
%                     I(:,:,j) = It(v(j)-r:v(j)+r,u(j)-r:u(j)+r);
%                 end
%                 sim(:,1) = sum(sum((I-Ip).^2));
%                 xt = u(sim==min(sim));
%                 yt = v(sim==min(sim));
%                 x(n) = round(xt(1))+xprev-rad;
%                 y(n) = round(yt(1))+yprev-rad; 
%             end
        end
    end
    
    %Plot cell tracks
    if i > 1
        delete(l3(:))
        delete(p(:))
        for j=1:size(cellPos,3)
            p(j) = plot(cellPos(:,1,j),cellPos(:,2,j),'Color',col(cellPos(end,4,j),:));
            l3(j) = line([x(j),cellPos(end,1,j)],[y(j),cellPos(end,2,j)],'Color',col(cellPos(end,4,j),:));
        end
    end

    %Scatter cells
    delete([s(:);t])
    [s,t] = ScatterCells(x,y,i-1,div);

    
    %Gets cell positions
    lock = zeros(numel(x),1);
    hide = 0;
    while 1
        clear midx dmin
        dmin = zeros;

        %Gets Input
        [xtemp,ytemp,btn] = ginput(1);

        if xtemp <= xres && xtemp >= 0 && ytemp <= yres && ytemp >= 0  
            for j=1:numel(x)
                dmin(j)=pdist([ytemp,xtemp;y(j,1) x(j,1)]);
            end
            dtemp = dmin;
            dmin(lock==1) = [];
            
            if isempty(dmin)
                midx = dtemp == min(dtemp);
            else
                midx = dtemp == min(dmin);
            end
            
            while sum(midx) > 1
                midx(find(midx,1)) = 0;
            end
            
            switch btn
                %Adds new cell - Left click
                case 1
                    n=n+1;
                    if i>1
                        cellPos(1:i-1,1,n) = cellPos(i-1,1,midx);
                        cellPos(1:i-1,2,n) = cellPos(i-1,2,midx);
                        cellPos(1:i-1,3,n) = cellPos(i-1,3,midx);
                        cellPos(i-1,4,n) = cellPos(i-1,4,midx);
                        div(n,1) = 0;
                    else
                        div(n,1) = 1;
                    end
                    
                    x(n,1) = round(xtemp);
                    y(n,1) = round(ytemp);

                %Repositions existing cell - Right click
                case 3
                    xp = x;
                    yp = y;
                    
                    try
                        x(midx) = round(xtemp);
                        y(midx) = round(ytemp);
                    catch
                    end
                
                %Undo - Ctrl Z    
                case 26
                    x=xp;
                    y=yp;
                    
                %Highlight selected track - B
                case 98
                    if bold
                        bold = 0;
                        delete([pb,lb])
                    else
                        bold = 1;
                        pb = plot(cellPos(:,1,midx),cellPos(:,2,midx),'Color','m','LineWidth',1.5);
                        lb = line([x(midx),cellPos(end,1,midx)],[y(midx),cellPos(end,2,midx)],'Color','m','LineWidth',1.5);
                    end
                    
                %Lock Cell Position - D
                case 100
                    lock(midx) = 1;
                    
                %Unlock Cell Position - F
                case 102
                    lock(midx) = 0;
                    
                %Hide Cell Tracks - H
                case 104
                    if hide == 0
                        hide = 1;
                        delete(l3(:))
                        delete(p(:))
                        delete(s(:))
                        try delete([pb,lb]); catch; end
                    end
                    
                %Show Cell Tracks - J
                case 106
                    if hide == 1
                        hide = 0;
                        for j=1:size(cellPos,3)
                            p(j) = plot(cellPos(:,1,j),cellPos(:,2,j),'Color',col(cellPos(end,4,j),:));
                            l3(j) = line([x(j),cellPos(end,1,j)],[y(j),cellPos(end,2,j)],'Color',col(cellPos(end,4,j),:));
                        end
                    end
                    
                case 115
                    if st
                        st=0;
                    else
                        st=1;
                    end
                %Deletes cell - Delete key
                case 127
                    try
                        x(midx) = [];
                        y(midx) = [];
                        cellPos(:,:,midx) = [];
                        delete([p(midx),l3(midx)])
                        try delete([pb,lb]); catch; end
                        n=n-1;
                    catch
                    end
            end
        end
        
        if btn==8 || btn==32
            hide = 0;
            bold = 0;
        end
        
        %Plot cell tracks
        if i > 1 && ~hide && ~bold
            delete(l3(:))
            delete(p(:))
            try delete([pb,lb]); catch; end
            for j=1:size(cellPos,3)
                p(j) = plot(cellPos(:,1,j),cellPos(:,2,j),'Color',col(cellPos(end,4,j),:));
                l3(j) = line([x(j),cellPos(end,1,j)],[y(j),cellPos(end,2,j)],'Color',col(cellPos(end,4,j),:));
            end
        end
        
        %Scatter cells
        delete([s(:);t])
        [s,t] = ScatterCells(x,y,i-1,div);
  
        %Goto Prev Image - Backspace
        if btn == 8 && i>1
            
            delete(l3(:))
            delete(p(:))
            try delete([pb,lb]); catch; end
            clear x y
            if i>2
                i = i-1;
                x(:,1) = cellPos(end,1,:);
                y(:,1) = cellPos(end,2,:);
                cellPos(end,:,:) = [];
                for j=1:size(cellPos,3)
                    if cellPos(end,4,j) == 0
                        cellPos(:,:,j) = [];
                        n=n-1;
                        div(j)=[];
                        x(j) = [];
                        y(j) = [];
                    end
                end
            else
                i=1;
                x(:,1) = cellPos(1,1,:);
                y(:,1) = cellPos(1,2,:);
            end
            break
        end
       
        
        %Goto Next Image - Space bar
        if btn == 32
            
            %Save Image
            for j=1:numel(x)
                yt(1:5) = m*x(j)+b+linspace(-1,3,5)*dx*rto;
                yt(6:10) = yt(1:5)+dx;
                yp = y(j);
                cellPos(i,1:2,j) = [x(j) y(j)];
                if (yp > yt(1) && yp < yt(6)) || (yp > yt(2) && yp < yt(7)) || (yp > yt(3) && yp < yt(8)) || (yp > yt(4) && yp < yt(9)) || (yp > yt(5) && yp < yt(10))
                    cellPos(i,3,j) = 1;
                else
                    cellPos(i,3,j) = 0;
                end
                if i > 1
                    cellPos(i,4,j) = cellPos(i-1,4,j);
                else
                    cellPos(i,4,j) = j;
                end
            end

            i = i+1;
            break
        end
    end
    
    %Save Cell Positions
    str = sprintf('%s/S%dData',fldrSave,fileNum);
    save(str,'cellPos')
end


%Scatter cell function
function [s,t] = ScatterCells(x,y,im,div)
    global m b dx rto
    s = gobjects;
    for i=1:numel(y)
        yt(1:5) = m*x(i)+b+linspace(-1,3,5)*dx*rto;
        yt(6:10) = yt(1:5)+dx;
        yp = y(i);
        if (yp > yt(1) && yp < yt(6)) || (yp > yt(2) && yp < yt(7)) || (yp > yt(3) && yp < yt(8)) || (yp > yt(4) && yp < yt(9)) || (yp > yt(5) && yp < yt(10))
            
            s(i,1) = scatter(x(i),y(i),10,'c','filled');
            if ~div(i)
                s(i,2) = scatter(x(i),y(i),3.5,'m ','filled'); 
            end
        else
            s(i,1) = scatter(x(i),y(i),10,'g','filled');
            if ~div(i)
                s(i,2) = scatter(x(i),y(i),3.5,'r','filled'); 
            end
        end
    end
    txt = sprintf('t=%0.2f h, n=%d',im/4,size(y,1));
    t = text(3,506,txt,'Color','white');
end

function Organize(old,numVids)
    global prefix fldrOpen
    for i=1:numVids

        %Create folder named 'S#' for each video position
        filename = sprintf('%s/S%d',fldrOpen,i);
        mkdir(filename)

        %Renames or moves image to correct folder
        for j=1:97
            oldName = sprintf('%s_s%d_t%d.tif',[fldrOpen,'/',old],i,j);
            newName = sprintf('%s_s%d_t%d.tif',[filename,'/',prefix],i,j);
            try
                movefile(oldName,newName)
            catch
                fprintf('%s Not Found\n',oldName)
            end
        end
    end
end


