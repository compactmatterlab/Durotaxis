clear global
clearvars
global dt kon F_thres kon_Act uInt pPax Fs v0

dt=1;
time = 0:dt:180;

%Active Stress Fibers (aSF) - binding of intergrin makes them "active"
numCells = 2;
numSF = 200;
maxSF = 75;


maxPa = 2500;
maxKeff = maxPa*0.1e-6; %N/m2 x m;

Fmotor = 2e-12;
motors = linspace(1,30,25); %

kon = 0.1;
F_thres = 1e-12;    %Force threshold (N)
kon_Act = 1;      %Actin-Paxillin on rate (s^-1)
uInt = 1.5;
pPax = 0.7;
v0 = 1e-6; %Velocity (m/s)
Keff =  [350*0.1e-6,1800*0.1e-6];linspace(0,maxKeff,numCells);
F = Fmotor.*motors;



SF_idx = zeros(1,numCells+1);
SF = struct('num',0,'bonds',0,'bound',0,'force',0,'size',0);
for i=1:numCells
    SF(i).num = numSF;
    SF(i).Keff = Keff(i)*ones(1,SF(i).num);
    SF_idx(i+1) = SF(i).num+SF_idx(i);
    
    SF(i).pos = [2*pi*(rand(1,SF(i).num));10/2*ones(1,SF(i).num);zeros(1,SF(i).num);zeros(1,SF(i).num)];
    SF(i).bonds = zeros(maxSF,SF(i).num,6);
    SF(i).bound = zeros(1,SF(i).num);
    SF(i).force = zeros(1,SF(i).num);
    SF(i).size = zeros(1,SF(i).num);
end
aSF = [SF.bonds];
K = [SF.Keff];

FASize = zeros(numel(time),numSF,numCells);
f = zeros(numel(time),numSF,numCells);
tavg = zeros(numSF,1);

G = numel(motors);
H = 20;
tb = zeros(numCells,G,H);
avgF = zeros(numCells,G,H);
avgF2 = zeros(numCells,G,H);
ton = zeros(numCells,G,H);

for g=1:G

   Fs = F(g);
    for h=1:H
       

       disp([g h])

        n1=1;
        for t=time
            aSF=Bind(K,aSF);
            for i=1:numCells
                SF(i).bonds = aSF(:,1+SF_idx(i):SF_idx(i+1),:);
                SF(i).bound = sum(SF(i).bonds(:,:,3))>0;
                SF(i).size = sum(SF(i).bonds(:,:,3));
                SF(i).force = sum(SF(i).bonds(:,:,1));

                FASize(n1,:,i) = SF(i).size;
                f(n1,:,i) = SF(i).force;
                Fm = SF(i).force;
                avgF2(i,g,h) = mean(mean(Fm(Fm>0)));
            end
            n1=n1+1;
        end

        for i=1:numCells
            for j=1:numSF
                tm = diff(find(FASize(:,j,i)==0));
                tm(tm==1) = [];
                tavg(j,1) = dt*nanmean(tm);

            end
            ton(i,g,h) = mean(sum(FASize(:,:,i)>0)./size(FASize,1));
            avgF(i,g,h) = mean(mean(f(:,:,i),2));
            tb(i,g,h) = nanmean(tavg);
        end
    end
end
save('F_vs_LT+F_Slip','avgF','tb','ton','Keff','motors','avgF2')

figure
c1 = [0,0.4470,0.7410];
c2 = [0.8500,0.3250,0.0980];
y1=mean(avgF2(1,:,:),3)*1e12;
y2=mean(avgF2(2,:,:),3)*1e12;
fancyPlot({motors*2,motors*2},{y1,y2},{'xlabel','Max SF Force (pN)'},...
    {'ylabel','SF Force (pN)'},{'ylim',[0,Inf]},{'xlim',[0,Inf]},...
    {'legend', 'Soft (0.35 kPa)','Stiff (1.8 kPa)'},{'spline'},...
    {'color',c1,c2},{'lineWidth',2},{'fontSize',18});

figure
y1=mean(tb(1,:,:),3);
y2=mean(tb(2,:,:),3);
fancyPlot({motors*2,motors*2},{y1,y2},{'xlabel','Max SF Force (pN)'},...
    {'ylabel','Bond Lifetime'},{'ylim',[0,Inf]},{'xlim',[0,Inf]},...
    {'spline'},{'color',c1,c2},{'lineWidth',2},{'fontSize',18});


% cl = get(groot,'defaultAxesColorOrder');
% figure
% set(gcf, 'Position', [20, 50, 800, 500])
% for i=1:numel(motors)
%     leg{i} = sprintf('%0.1f pN',motors(i)*2);
%     y1=nanmean(avgF2(:,i,:),3)*1e12;
%     y1(isnan(y1)) = 0;
% 
%     pl(i) = fancyPlot({Keff./1e-4},{y1},{'xlabel','Elasticity (kPa)'},...
%         {'ylabel','SF Force (pN)'},{'ylim',[0,100]},{'xlim',[0,Inf]},...
%          {'spline'},{'color',cl(mod(i-1,7)+1,:)},{'lineWidth',2},{'fontSize',18}); 
% end
% l=legend(pl(1:end),leg{1:end},'Location','southeast');
% l.EdgeColor = [1,1,1];
% l.FontSize = 18;
% set(l,'EdgeColor','none');
% set(l,'color','none');
% 
% figure
% set(gcf, 'Position', [20, 50, 800, 500])
% for i=1:numel(motors)
%     leg{i} = sprintf('%0.1f pN',motors(i)*2);
%     y1=mean(tb(:,i,:),3);
% 
%     pl(i) = fancyPlot({Keff./1e-4},{y1},{'xlabel','Elasticity (kPa)'},...
%         {'ylabel','Average Bond Lifetime (s)'},{'ylim',[0,Inf]},{'xlim',[0,2]},...
%         {'spline'},{'color',cl(mod(i-1,7)+1,:)},{'lineWidth',2},{'fontSize',18});
% end
% l=legend(pl(1:end),leg{1:end},'Location','southeast');
% l.EdgeColor = [1,1,1];
% l.FontSize = 18;
% set(l,'EdgeColor','none');
% set(l,'color','none');


function SFbonds = Bind(Keff,aSF)
global dt kon F_thres kon_Act uInt pPax


Fsf = aSF(:,:,1);       %Force on Stress Fiber
tsf = aSF(:,:,2);       %Time force has been genereated by SF
Nt = aSF(:,:,3);        %Number of available binding sites
Nb = aSF(:,:,4);        %Number of attached bonds
Pax = aSF(:,:,5);       %Stress sensor protein 
SFprev = aSF(:,:,6);    
maxSF = size(Fsf,1);

%Generate # of binding sites for unbound SFs
idx = ~sum(Nt);
if sum(idx)
    Nt(1,idx) = fpoissrnd(uInt,sum(idx),1);
    Pax(1,idx) = fbinornd(Nt(1,idx),pPax);
end


%Test bonds
dN = fbinornd(Nt-Nb,1-exp(-dt.*kon))-fbinornd(Nb,1-exp(-dt.*kOff_slip(Fsf./Nb)));
Nb = Nb + dN;


%Reopen paxillins of parent SFs and delete unbound SFs
idx = SFprev>0 & ~Nb;
rePax = SFprev(idx);
for j=1:numel(rePax)
    Pax(rePax(j)) = Pax(rePax(j))+1;
end
tsf(~Nb) = 0;
Nt(~Nb) = 0;
Pax(~Nb) = 0;
SFprev(~Nb) = 0;

%Generate force
Fsf = Force(Keff,tsf);
tsf(Nb>0) = tsf(Nb>0) + dt;

%Add new aSF to available Paxillins
psf = find(Fsf./Nb>F_thres);
new = fbinornd(Pax(psf),1-exp(-tsf(psf).*kon_Act));
Pax(psf) = Pax(psf) - new;
col = zeros(1,sum(new));
if ~isempty(col)
    n=1;
    for j=1:numel(new)
        for k=1:new(j)
            col(n) = ceil(psf(j)/maxSF);
            n=n+1;
        end
    end
    
    idx = ~Nt;
    row = zeros(size(col));
    for j=1:numel(col)
        try
            row(j) = find(idx(:,col(j)),1);
            idx(row(j),col(j)) = 0;
        catch
            row(j) = row(j-1);
        end
    end
    add = sub2ind(size(Nt),row,col)';
    Nt(add) = fpoissrnd(uInt,numel(add),1);
    Pax(add) = fbinornd(Nt(add),pPax);
    prev = zeros(1,sum(new));
    n=1;
    for j=1:numel(new)
        for k=1:new(j)
            prev(n) = psf(j);
            n=n+1;
        end
    end
    SFprev(add) = prev;
end

SFbonds(:,:,1) = Fsf;
SFbonds(:,:,2) = tsf;
SFbonds(:,:,3) = Nt;
SFbonds(:,:,4) = Nb;
SFbonds(:,:,5) = Pax;
SFbonds(:,:,6) = SFprev;
end


function k = kOff(f)
f(isnan(f)) = 0;
kb = 1.3806e-23;
T = 310;
fx = 5.3800e-12;
xi = kb*T/fx;

%a5B1 - FN (Ca+Mg)
a = 3.309;
b = 0.0003942;
c = 58.19;

% % %a5B1 - FN (Mg Only)
% a = 1.5774;
% b = 0.0003;
% c = 14.3570;

% %Truncated a5B1 - FN (Ca+Mg)
% a = 3.6225;
% b = 3.7605e-04;
% c = 3.0726;

%This equation I made to fit the data from Kong 2009, JCB (Inversed to go from <t> to koff)
k = (a*exp(-f*xi./(kb*T))+(b*exp(f*xi./(kb*T))+c*exp(-f*xi./(kb*T))).^-1).^-1;


end

function k = kOff_slip(f)

k0 = 0.1;
Fb = 2e-12;

%Bell equation with constants from Bangasser/Odde 2013 Cell Mol Bioeng.
k = k0*exp(f./Fb);
end


function f=Force(Keff,t)
global Fs v0

f=Fs.*(1-exp(-(v0*Keff.*t)./Fs));
end
