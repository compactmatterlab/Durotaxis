clear all
scale = 7.5879; %pixels/um
E = 'Soft';
A = 'WA';

E_all = {'Stiff' 'Soft' 'Stiff' 'Soft'};
A_all = {'SA' 'SA' 'WA' 'WA'};
for h=1:4
    E = E_all{h};
    A = A_all{h};
    n=1;
for j=1:3
for i=1:10
    image = [A int2str(j) ' ' E '/', 'MDA-MB-231_' A int2str(j) '_' E '_Mar32020-', int2str(i)];
    
    try
    %% Create mask from actin
    TXRD = imread([image, '_TXRD.tif']);
    I = TXRD;
    T = adaptthresh(I,0.2,'Statistic','gaussian');
    M = imbinarize(I,T);
    M = imdilate(M,strel('disk',1));
    CC = bwconncomp(M);
    S = regionprops(CC, 'Area');
    area = round(max([S.Area]))-1;
    M = bwareaopen(M, area);
    M = imclose(M,strel('disk',8,4));
    M = imfill(M,'holes');
%     subplot(2,5,i)
%     imshow(M)
  
%     %% Create mask from cell mask
%     DPRD = imread([image, '_DPRD.tif']);
%     I = DPRD;
%     T = adaptthresh(I,0.3,'Statistic','gaussian');
%     M1 = imbinarize(I,T);
%     CC = bwconncomp(M1);
%     S1 = regionprops(CC, 'Area');
%     area = round(max([S.Area]))-1;
%     M1 = bwareaopen(M1, area);
%     M1 = imclose(M1,strel('disk',8,4));
%     M1 = imfill(M1,'holes');
%     
% %     subplot(2,5,i)
%     imshow(M1)
    
    %% Get Focal Adhesions
    if h~=2
     clear I
    end
    FITC = imread([image, '_FITC.tif']);
    I(:,:,1) = double(rescale(FITC));
    for k=1:3
        iTemp = I(:,:,k);
        T = prctile(iTemp(:),99);
        M2 = imbinarize(I(:,:,k),T);
        M2 = M&M2;
        M2 = bwareaopen(M2, 3); %Min area cutoff (50 nm2)
        CC = bwconncomp(M2);
        S2 = regionprops(CC, 'Area');
        try
            I(:,:,k+1) = I(:,:,k) - bwareafilt(M2,sum([S2.Area]>575)); %Max area cutoff (10um2)
        catch
            I(:,:,k+1) = I(:,:,k);
        end
        iTemp = I(:,:,k+1);
        T = prctile(iTemp(:),99.5);
        try
            M2 = imbinarize(M2,T);
        end
        M2 = M&M2;
        M2 = bwareaopen(M2, 3); %Min area cutoff (50 nm2)
        CC = bwconncomp(M2);
        S2 = regionprops(CC, 'Area');
        try
            M2 = M2 - bwareafilt(M2,sum([S2.Area]>575)); %Max area cutoff (10um2)
        catch
        end
    end
    
    %Get FA areas within mask
    CC = bwconncomp(M2);
    S3 = regionprops(CC, 'Area');
    FA = [S3.Area]./scale^2;
    avgFA(n,h) = rmoutliers(mean(FA));
    numFA(n,h) = numel(FA);
    
    CC = bwconncomp(M);
    S4 = regionprops(CC, 'Area');
    cellArea(n,h) = [S4.Area]./scale^2;
%     subplot(5,6,n)
%     imshow(M2)
    try
        O(:,:,:,n,h) = labeloverlay(FITC,M2);
    end
    
    catch e
        fprintf(2,'%s\n',e.message);
        avgFA(n,h) = NaN;
        numFA(n,h) = NaN;
        cellArea(n,h) = NaN;
    end
    n=n+1;
end
end
end
%implay(O,8)

col=[0,0.447,0.741;
    0.30196,0.7451,0.93333;
    0.85098,0.32549,0.098039;
    0.87059,0.4902,0;
    0.8,0.8,0.8;
    0.74902,0,0.74902;
    0.46667,0.67451,0.18824];

subplot(1,2,1)
[h1,p1]=beeswarm({avgFA(:,3),'WA Stiff',col(1,:)},{avgFA(:,4),'WA Soft',col(2,:)},...
    {avgFA(:,1),'SA Stiff',col(3,:)},{avgFA(:,2),'SA Soft',col(4,:)},...
    {'Title',' '},{'ylabel','FA Area (\mum^2)'},{'stats'},{'Outliers'},{'FontSize',16});

subplot(1,2,2)
numFA(numFA==0)=NaN;
[h2,p2]=beeswarm({numFA(:,3)./cellArea(:,3),'WA Stiff',col(1,:)},{numFA(:,4)./cellArea(:,4),'WA Soft',col(2,:)},...
    {numFA(:,1)./cellArea(:,1),'SA Stiff',col(3,:)},{numFA(:,2)./cellArea(:,2),'SA Soft',col(4,:)},...
    {'Title',' '},{'ylabel','# FA / Cell Area (\mum^-^2)'},{'stats'},{'Outliers'},{'FontSize',16});

figure
imshow(O(:,:,:,13,1))
figure
imshow(O(:,:,:,12,2))
figure
imshow(O(:,:,:,14,3))
figure
imshow(O(:,:,:,11,4))