function DrawCells(cellPos,t,type,step,aSF,pConsts,mParams,video,SF_idx)
global lim limx limy lnew

dt = mParams.dt;
numCells = mParams.numCells;

% if t==1 && type~=2
%     set(gcf, 'Position',[100 0   840   792.2])
% end
if mod(t-1,step*60)==0 && video~=0
    
    wSoft = 170;
    wStiff = 115;

    clf
    %Get direction of polarity
    dir = zeros(2,numCells);
    dir(1:numel([aSF.dir])) = [aSF.dir];
    pTheta = mod(cart2pol(dir(1,:),dir(2,:)),2*pi);
    for c=1:numCells
        idx = SF_idx(c)+1:SF_idx(c+1);
        SFpos = (aSF.pos(:,idx))';
        SFbonds = aSF.bonds(:,idx,:);
        force = aSF.force(idx)';
        faSize = aSF.size(idx)';
        cellDiameter = aSF.cDiameter(c);
        
        numSF = sum(SFbonds(:,:,3))';
        
        iLead = aSF.lead;%(mod(SFpos(:,1)-pTheta(c),2*pi))<pi/2 | (mod(SFpos(:,1)-pTheta(c),2*pi))>=3*pi/2;

        [~,I] = sort(SFpos(:,1)-min(SFpos(:,1)));
        dc_sFibers = SFpos(I,:);
        dc_faSize = faSize(I);
        dc_iLead = iLead(I);
        dc_force = force(I,:).*SFpos(I,3:4)./vecnorm(SFpos(I,3:4));
        dc_numSF = numSF(I);
        dc_cellPos = cellPos(:,1:t,c)';
        
        sFibersL = dc_sFibers(dc_iLead,:,:);
        if max(sFibersL(:,1))-min(sFibersL(:,1))>pi
            [~,I] = sort(mod(sFibersL(:,1)+pi,2*pi));
            sFibersL = sFibersL(I,:);
        end
        sFibersT = dc_sFibers(~dc_iLead,:);
        if max(sFibersT(:,1))-min(sFibersT(:,1))>pi
            [~,I] = sort(mod(sFibersT(:,1)+pi,2*pi));
            sFibersT = sFibersT(I,:);
        end
        faSizeL = dc_faSize(dc_iLead);
        faSizeT = dc_faSize(~dc_iLead);
        forceL = dc_force(dc_iLead,:);
        forceT = dc_force(~dc_iLead,:);
        numSFL = dc_numSF(dc_iLead);
        numSFT = dc_numSF(~dc_iLead);       
        
        switch type
            case 0
                set(gcf, 'Position',[100 0   840   792.2])
                X1 = [wSoft/2 wSoft/2+wStiff wSoft/2+wStiff wSoft/2];
                X2 = X1+(wSoft+wStiff);
                patch(-X1,[-600 -600 600 600],[0.01,0.01,0.01], 'FaceAlpha',0.1, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(X1,[-600 -600 600 600],[0.01,0.01,0.01], 'FaceAlpha',0.1, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(X2,[-600 -600 600 600],[0.01,0.01,0.01], 'FaceAlpha',0.1, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(-X2,[-600 -600 600 600],[0.01,0.01,0.01], 'FaceAlpha',0.1, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                title(sprintf('F=%0.1f pN',nanmean(pConsts.Fmax*1e12)))
                text(65*1.02,-20*.95,sprintf('t=%0.2f hr',dt*(t-1)/3600))
                Cell
%                 Trajectory
%                 title(sprintf('F=%d pN',pConsts.Fmax*1e12))
%                 text(-lim*.95,-lim*.95,sprintf('t=%0.2f hr',dt*(t-1)/3600))
%                 axes('Position',[.7 .7 .2 .2])
%                 Cell
                
             case 1
                set(gcf, 'Position',[100 0   840   792.2])
                Cell
                title(sprintf('F=%d pN',pConsts.Fmax*1e12))
                text((-19+limx),(-19+limy),sprintf('t=%0.2f hr',dt*(t-1)/3600))
                axes('Position',[.7 .7 .2 .2])
                Trajectory
            case 2
                if c==1; fprintf('t=%0.2f hr\n',dt*(t-1)/3600); end
                
            otherwise
                set(gcf, 'Position',[100 0   840   792.2])
                X1 = [wSoft/2 wSoft/2+wStiff wSoft/2+wStiff wSoft/2];
                X2 = X1+(wSoft+wStiff);
                patch(-X1,[-600 -600 600 600],[0.1,0.1,0.1], 'FaceAlpha',0.01, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(X1,[-600 -600 600 600],[0.1,0.1,0.1], 'FaceAlpha',0.01, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(X2,[-600 -600 600 600],[0.1,0.1,0.1], 'FaceAlpha',0.01, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                patch(-X2,[-600 -600 600 600],[0.1,0.1,0.1], 'FaceAlpha',0.01, 'EdgeColor', [0.2,0.2,0.2], 'EdgeAlpha', 1, 'LineWidth', 0.5);
                Trajectory
                title(sprintf('F=%0.1f pN',mean(pConsts.Fmax*1e12)))
                if c==1
                    text(-lim*.95,-lim*.95,sprintf('t=%0.2f hr',dt*(t-1)/3600))
                end
                
        end
    end
    if type~=2
        drawnow
        try
            frame = getframe(gcf);
            writeVideo(video,frame);
        catch
        end
    end
end
    function Trajectory
        if t<=1
            lnew = 100;
        end
        
        ltemp = max([ceil(max(abs(dc_cellPos(t,1)))/100)*100,ceil(max(abs(dc_cellPos(t,2)))/100)*100]);
        if ltemp > lnew
            lnew = ltemp;
        end
        
        box on
        hold on
        plot(dc_cellPos(1:step*60:end,1),dc_cellPos(1:step*60:end,2))
        try
            [pL(1,1),pL(1,2)] = pol2cart(sFibersL(1,1),cellDiameter/2);
            [pL(2,1),pL(2,2)] = pol2cart(sFibersL(end,1),cellDiameter/2);
            posL = [pL(1,:);sFibersL(:,3:4);pL(2,:)]+dc_cellPos(t,:);
              %posL = [sFibersL(:,3:4);sFibersT(1,3:4);]+dc_cellPos(t,:);
        catch
            posL = dc_cellPos(t,:);
        end
        try
            [pT(1,1),pT(1,2)] = pol2cart(sFibersT(1,1),cellDiameter/2);
            [pT(2,1),pT(2,2)] = pol2cart(sFibersT(end,1),cellDiameter/2);
            posT = [pT(1,:);sFibersT(:,3:4);pT(2,:)]+dc_cellPos(t,:);
              %posT = [sFibersT(:,3:4);sFibersL(1,3:4)]+dc_cellPos(t,:);
        catch
            posT = dc_cellPos(t,:);
        end
        
        plot(posL(:,1),posL(:,2),'color','k');
        plot(posT(:,1),posT(:,2),'color','k');
        circle(dc_cellPos(t,1),dc_cellPos(t,2),cellDiameter/2,'color','k','filled');
        
        lim = 300;
        xlim([-lim,lim])
        ylim([-lim,lim])
%         if t<=1
%             lim = 100;
%             xlim([-lim,lim])
%             ylim([-lim,lim])
%         else
%             lim = lnew;
%             xlim([-lim,lim])
%             ylim([-lim,lim])
%         end
    end


    function Cell
        box on
        hold on
        
        posL = sFibersL(:,3:4)+dc_cellPos(t,:);
        posT = sFibersT(:,3:4)+dc_cellPos(t,:);
        if ~isempty(posL)
            s1 = scatter(posL(faSizeL>0,1),posL(faSizeL>0,2),45*(faSizeL(faSizeL>0)+0.001),[0.03,0.03,0.03],'filled');
        end
        if ~isempty(posT)
            s2 = scatter(posT(faSizeT>0,1),posT(faSizeT>0,2),45*(faSizeT(faSizeT>0)+0.001),[0.2,0.2,0.2],'filled');
        end
        for i=1:size(posL,1)
            line([dc_cellPos(t,1),posL(i,1)],[dc_cellPos(t,2),posL(i,2)],'color',s1.CData,'linewidth',0.01+0.02*numSFL(i));
        end
        for i=1:size(posT,1)
            line([dc_cellPos(t,1),posT(i,1)],[dc_cellPos(t,2),posT(i,2)],'color',s2.CData,'linewidth',0.01+0.02*numSFT(i));
        end
        
       
        %plot(dc_cellPos(:,1),dc_cellPos(:,2))
        clear posL posT
        try
            [pL(1,1),pL(1,2)] = pol2cart(sFibersL(1,1),cellDiameter/2);
            [pL(2,1),pL(2,2)] = pol2cart(sFibersL(end,1),cellDiameter/2);
            posL = [pL(1,:);sFibersL(:,3:4);pL(2,:)]+dc_cellPos(t,:);
        catch
            posL = dc_cellPos(t,:);
        end
        
        
        try
            [pT(1,1),pT(1,2)] = pol2cart(sFibersT(1,1),cellDiameter/2);
            [pT(2,1),pT(2,2)] = pol2cart(sFibersT(end,1),cellDiameter/2);
            posT = [pT(1,:);sFibersT(:,3:4);pT(2,:)]+dc_cellPos(t,:);
        catch
            posT = dc_cellPos(t,:);
        end
        
        plot(posL(:,1),posL(:,2),'color','k');
        plot(posT(:,1),posT(:,2),'color','k');
        
        circle(dc_cellPos(t,1),dc_cellPos(t,2),cellDiameter/2,'color','k','filled');
        F = 3*(sum(forceL,1)+sum(forceT,1));
        quiver(dc_cellPos(t,1),dc_cellPos(t,2),F(1),F(2),1e11,'color','c','linewidth',3,'MaxHeadSize',0.5)
        quiver(dc_cellPos(t,1),dc_cellPos(t,2),F(1),0,1e11,'color',[128/255,1,0],'linewidth',3,'MaxHeadSize',0.5)
        quiver(dc_cellPos(t,1),dc_cellPos(t,2),0,F(2),1e11,'color','r','linewidth',3,'MaxHeadSize',0.5)
        
        
        limx = dc_cellPos(t,1);
        limy = dc_cellPos(t,2);
        xlim([-20+limx,20+limx])
        ylim([-20+limy,20+limy])
    end


end