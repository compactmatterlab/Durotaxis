%% 2D Cell Migration Model for Simulating Durotaxis
%code written by Ben Yeoman
%created 08/22/2019
%edited 05/07/2020
clear global
clearvars
close all

%**************************************************************************
%% Model parameters
mParams.runTime = 3600*24;
mParams.dt = 1;
mParams.numCells = 10;
mParams.samples = 10;
mParams.seed = sum(clock);
saveAll = 1;
video = 0;
rng(mParams.seed)

Vars(1,:)=linspace(1e-18,70e-12,mParams.samples);
Vars(2,:)=linspace(0,3,mParams.samples);
Vars(3,:)=linspace(0,1,mParams.samples);

%**************************************************************************
%% Physical constants
pConsts.Fmotor = 2e-12;                             %Stall force for a single myosin II motor (N)
pConsts.motors = 45/2;                              %Average number of myosin II motors per F-actin-
pConsts.Fmax = pConsts.Fmotor*pConsts.motors;  %%%  %Max Stress Fiber force (N) 
pConsts.friction = 2e-5;                            %Bond friction (m * kg m-1 s-1)
pConsts.tRtrc = 10;                                 %Average time to retraction (s)
pConsts.F_thres = 1e-12;                            %Force sensor threshold (N)
pConsts.kon = 0.1;                                  %Integrin-SF assembly rate (s^-1)
pConsts.kon_Act = 1;                                %Actin-Paxillin on rate (s^-1)
pConsts.uInt = 1;                              %%%  %Average Integrins/F-actin
pConsts.pPax = 0.7;                            %%%  %Probability of Paxillin bound to Integrin
pConsts.vRet = 0.5;                                 %Average SF retraction velocity (um/s)
pConsts.vExtL = 0.2;                                %Average leading edge SF growth velocity (um/s)
pConsts.vExtT = 0.1;                                %Average trailing edge SF growth velocity (um/s)
pConsts.Estiff = 30000;                              %Elasticity of stiff gel (Pa)
pConsts.Esoft = 10000;                                %Elasticity of soft gel (Pa)
pConsts.Kstiff = pConsts.Estiff*0.1e-6;             %Stiffness of stiff gel (N/m)
pConsts.Ksoft = pConsts.Esoft*0.1e-6;               %Stiffness of soft gel (N/m)
pConsts.Gradient = 10;                              %Length of gradient between soft and stiff (um) 
pConsts.maxSF = 100;                                %Limit for numer of stress fiber/seed
pConsts.avgSF = 50;                                 %Avg number of initiation sites/cell
pConsts.v0 = 1e-6;                                  %Maximum (unloaded) myosin II sliding velocity on F-actin (m/s)

%**************************************************************************
%% Saving and video setup 
if saveAll
    
    %Saving
    saveAs = sprintf('Migration2D_n=%d_Fs=%dpN_1',mParams.numCells*mParams.samples,pConsts.Fmax*1e12); %#ok<*UNRCH>
    i=2;
    while isfile(['Results\',saveAs, '.mat'])
        saveAs = [saveAs(1:end-numel(num2str(i-1))),num2str(i)];
        i=i+1;
    end

    %Movie Maker
    if video
        vidName = ['Movies\', saveAs];
        video = VideoWriter(vidName);
        video.FrameRate = 8;
        open(video);
    end
end

%**************************************************************************
%% Sliced structures for parallel computing

%Model Parameters
runTime = mParams.runTime; 
dt = mParams.dt;
numCells = mParams.numCells;
samples = mParams.samples;

%Acto-Myosin Constants
Fmotor = pConsts.Fmotor;
motors = pConsts.motors;

%Bind Constants
kon = pConsts.kon;
F_thres = pConsts.F_thres;
kon_Act = pConsts.kon_Act;
uInt = pConsts.uInt;
pPax = pConsts.pPax;
Kstiff = pConsts.Kstiff;
Ksoft = pConsts.Ksoft;
maxSF = pConsts.maxSF;
avgSF = pConsts.avgSF; 

%Extend Constants
vRet = pConsts.vRet;
vExtL = pConsts.vExtL;
vExtT = pConsts.vExtT;
tRtrc = pConsts.tRtrc;

%Move Constants
friction = pConsts.friction; 

%**************************************************************************
%% Prealloctions for output parameters
cellPos = zeros(2,runTime/dt+2,numCells,samples);
Traction = zeros(numCells,runTime/dt+2,samples);
vel = zeros(numCells,samples);
dis = zeros(numCells,samples);

%**************************************************************************
%%
tic
for n=1:samples          
    %% Set experimental variables
    %var = [Vars(3,n),3];
    %Fmax = 0;
    
    %% Set initial conditions
    cPos = zeros(2,runTime/dt+2,numCells);
    tractF = zeros(numCells,runTime/dt+2);
    if numCells==1
        cPos(:,1,:) = [0,0];
    else
        cPos(1,1,:) = linspace(-285,285,numCells);%rand(1,1,numCells)*500-500/2;%
        cPos(2,1,:) = rand(1,1,numCells)*500-500/2;
    end
    
    SF = struct('num',0,'cDiameter',0,'pos',0,'bonds',0,'bound',0,'force',0,'size',0,'tRet',0,'dir',0,'pTheta',0);
    SF_idx = zeros(1,numCells+1);
    for i=1:numCells
        
        SF(i).num = fpoissrnd(avgSF);
        SF(i).cDiameter = fnormrnd(10,0.1,1,SF(i).num);
        
        SF(i).pos = [2*pi*(rand(1,SF(i).num));10/2*ones(1,SF(i).num);zeros(1,SF(i).num);zeros(1,SF(i).num)];
        SF(i).bonds = zeros(maxSF,SF(i).num,6);
        SF(i).bound = false(1,SF(i).num);
        SF(i).force = zeros(1,SF(i).num);
        SF(i).lead = zeros(1,SF(i).num);
        SF(i).size = zeros(1,SF(i).num);
        SF(i).tRet = -tRtrc.*log(rand(1,SF(i).num));

        SF(i).dir = rand(2,1); 
        SF(i).dir = SF(i).dir./norm(SF(i).dir); 
        
        SF_idx(i+1) = SF(i).num+SF_idx(i);
    end
    aSF.cDiameter = [SF.cDiameter];
    aSF.pos = [SF.pos];
    aSF.bonds = [SF.bonds];
    aSF.bound = [SF.bound];
    aSF.force = [SF.force];
    aSF.lead = [SF.lead];
    aSF.size = [SF.size];
    aSF.tRet = [SF.tRet];
    aSF.dir = [SF.dir];
    
    Fmax = fpoissrnd(motors,1,size(aSF.size,2)).*Fmotor;
%     nIdx = ceil(numel(SF_idx)/2);
%     Fmax(1:SF_idx(nIdx)) = fpoissrnd(motors,1,SF_idx(nIdx)).*Fmotor;
%     Fmax(SF_idx(nIdx)+1:SF_idx(end)) = fpoissrnd(30/2,1, numel(SF_idx(nIdx)+1:SF_idx(end))).*Fmotor;
    
    %**********************************************************************
    %% Main loop
    disp(n)
    t = 1;
    for rn=0:dt:runTime  
        aSF = Extend(aSF,pConsts,mParams,SF_idx);
        aSF = Bind(cPos(:,t,:),aSF,pConsts,mParams,SF_idx,Fmax);
        [cPos(:,t+1,:),aSF,tractF(:,t)] = Move(cPos(:,t,:),aSF,pConsts,mParams,SF_idx);
        DrawCells(cPos,t,1-numCells,15/dt,aSF,pConsts,mParams,video,SF_idx)
        if mod(t-1,15*60/dt)==0 && video==0
            fprintf('t=%0.2f \n',dt*t/3600)
        end
        t=t+1;
    end
    cellPos(:,:,:,n) = cPos; 
    Traction(:,:,n) = tractF;
    %**********************************************************************
end
dToc;

%**************************************************************************
%% Trajectory analysis
for n=1:samples
    for i=1:numCells
        [vel(i,n),dis(i,n)] = speed(cellPos(:,linspace(1,runTime/dt+1,(runTime)/900+1),i,n),runTime/3600);
        %[D(i,n),alp(i,n),R2(i,n)] = MSD(cPos(:,:,i),mParams.runTime,mParams.dt,1);
        %[Lp(i,n),R1(i,n)] = persistence(cPos(:,:,i),0);
    end
end
    
%**************************************************************************
%% Save data
if saveAll
    save(['Results\', saveAs], 'vel','dis','cellPos','Traction','mParams','pConsts')
    if video~=0
        close(video)
    end
end

%**************************************************************************
%% Stress fiber extension (or retraction)
function aSF = Extend(aSF,pConsts,mParams,SF_idx)
      
dt = mParams.dt;
numCells = mParams.numCells;

nu = 6.9130e-04;            %Viscosity of water at 37 C
temp = 310.15;              %Temperature of cell (K)
kb = 1.3806e-23;            %Boltzmann's constant
dAct = 7e-9;                %Diameter of actin filament (m)

vRet = pConsts.vRet;
vExtL = pConsts.vExtL;
vExtT = pConsts.vExtT;
tRtrc = pConsts.tRtrc;

%Get direction of polarity
dir = aSF.dir;
pTheta = zeros(1,SF_idx(end));
for i=1:numCells
    pTheta(SF_idx(i)+1:SF_idx(i+1)) = mod(cart2pol(dir(1,i),dir(2,i)),2*pi);
end

%Set variables
SFpos = aSF.pos;
tRet = aSF.tRet;
Ibound = aSF.bound;
cellDiameter = aSF.cDiameter;

%Determine leading and trailing edge
iLead = ((mod(SFpos(1,:)-pTheta,2*pi))<pi/2 | (mod(SFpos(1,:)-pTheta,2*pi))>=3*pi/2);

%Add rotational drift to stress fibers
idxL = iLead&~Ibound;
idxT = ~iLead&~Ibound;
lL = SFpos(2,idxL);
lT = SFpos(2,idxT);
DsfL = (3*kb*temp*reallog(lL./dAct))./(pi*nu*lL.^3);
DsfT = (3*kb*temp*reallog(lT./dAct))./(pi*nu*lT.^3);
SFpos(1,idxL) = SFpos(1,idxL) + fnormrnd(zeros(1,sum(idxL)),ones(1,sum(idxL)).*DsfL*dt);
SFpos(1,idxT) = SFpos(1,idxT) + fnormrnd(zeros(1,sum(idxT)),ones(1,sum(idxT)).*DsfT*dt);
    
%Extend or shrink stress fiber 
tRet(~Ibound) = tRet(~Ibound)-dt;
tRet(Ibound) = -tRtrc.*log(rand(1,numel(tRet(Ibound))));
    
%Retract 
idx = tRet<=0;
SFpos(2,idx) = SFpos(2,idx) - fnormrnd(dt*vRet*erf(SFpos(2,idx)-cellDiameter(idx)/2),0.1);

%Extend Leading
idx =~Ibound&iLead;
SFpos(2,idx) = SFpos(2,idx) + fnormrnd(dt*vExtL,0.01);

%Extend trailing
idx =~Ibound&~iLead;
SFpos(2,idx) = SFpos(2,idx) + fnormrnd(dt*vExtT,0.01);

%Get cartesian coordinates of stress fibers
[SFpos(3,:),SFpos(4,:)] = pol2cart(SFpos(1,:),SFpos(2,:));
SFpos(3:4,:) = round(SFpos(3:4,:),10);
    
%Update SF structures
aSF.pos = SFpos;
aSF.tRet = tRet;
aSF.lead = iLead;
end


%**************************************************************************
%% Focal Adhesion Formation
function aSF = Bind(curPos,aSF,pConsts,mParams,SF_idx,Fmax)

dt = mParams.dt;
numCells = mParams.numCells;

uInt = pConsts.uInt;
pPax = pConsts.pPax;

kon = pConsts.kon;
F_thres = pConsts.F_thres;
kon_Act = pConsts.kon_Act;
Kstiff = pConsts.Kstiff;
Ksoft = pConsts.Ksoft;
grad = pConsts.Gradient;
v0 = pConsts.v0;
maxSF = pConsts.maxSF;

%Get position for each stress fiber
pos = zeros(1,SF_idx(end));
for i=1:numCells
    pos(SF_idx(i)+1:SF_idx(i+1)) = aSF.pos(3,SF_idx(i)+1:SF_idx(i+1))+curPos(1,1,i);
end

%Concatanate all SFs for faster manipulation
Fsf = aSF.bonds(:,:,1);
tsf = aSF.bonds(:,:,2);
Nt = aSF.bonds(:,:,3);
Nb = aSF.bonds(:,:,4);
Pax = aSF.bonds(:,:,5);
SFprev = aSF.bonds(:,:,6);
    
%Generate # of binding sites for unbound SFs only
unbound = ~sum(Nt);
if sum(unbound) 
    Nt(1,unbound) = fpoissrnd(uInt,sum(unbound),1);
    Pax(1,unbound) = fbinornd(Nt(1,unbound),pPax);
end

%Test bonds
dN = fbinornd(Nt-Nb,1-exp(-dt.*kon))-fbinornd(Nb,1-exp(-dt.*kOff(Fsf./Nb)));
Nb = Nb + dN;

%Reopen paxillins of parent SFs and delete unbound SFs 
idxDelete = SFprev>0 & ~Nb;
rePax = SFprev(idxDelete);
for j=1:numel(rePax)
    Pax(rePax(j)) = Pax(rePax(j))+1;
end
tsf(~Nb) = 0;
Nt(~Nb) = 0;
Pax(~Nb) = 0;
SFprev(~Nb) = 0;

%Generate force based on position
Fsf = Force(Fmax,tsf,pos,Kstiff,Ksoft,grad,v0);
tsf(Nb>0) = tsf(Nb>0) + dt;

%Add new SF to available Paxillins
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

    idxN = ~Nt;
    row = zeros(size(col)); 
    for j=1:numel(col)
       try
           row(j) = find(idxN(:,col(j)),1);
           idxN(row(j),col(j)) = 0;
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

%Convert back to aSF
aSF.bonds(:,:,1) = Fsf;
aSF.bonds(:,:,2) = tsf;
aSF.bonds(:,:,3) = Nt;
aSF.bonds(:,:,4) = Nb;
aSF.bonds(:,:,5) = Pax;
aSF.bonds(:,:,6) = SFprev;

%Calculate bound, size, and force
aSF.bound = sum(Nt)>0;
aSF.size = sum(Nb);
aSF.force = sum(Fsf);

end

%**************************************************************************
%% Cell Movement
function [newPos,aSF,tracF] = Move(oldPos,aSF,pConsts,mParams,idx)

dt = mParams.dt;
numCells = mParams.numCells;

newPos = zeros(size(oldPos));
%tracF = zeros(size(oldPos));
tracF = zeros(numCells,1);
dir = zeros(2,1);
for i=1:numCells
    SFpos = aSF.pos(:,idx(i)+1:idx(i+1));
    Ibound = aSF.bound(idx(i)+1:idx(i+1));
    force = aSF.force(idx(i)+1:idx(i+1));
    %lead = aSF.lead(idx(i)+1:idx(i+1));
    numI = aSF.size(idx(i)+1:idx(i+1));
    
    friction = pConsts.friction*sum(numI); 
    
    %Get traction forces at leading and trailing ends
    %tracF(:,1,i) = [sum(force(lead))./sum(Ibound(lead));sum(force(~lead))./sum(Ibound(~lead))];
    tracF(i,1) = sum(force);
    
    %Calculate distance to move cell
    force = aSF.force(idx(i)+1:idx(i+1)).*SFpos(3:4,:)./vecnorm(SFpos(3:4,:));
    force(isnan(force)) = 0;
    F = norm(sum(force,2));

    if F
        dir = sum(force,2)./F;
        v = F/friction;
        d = v*dt*1e6;
        newPos(:,1,i) = d*dir+oldPos(:,1,i);
    else
        d=0;
        newPos(:,1,i) = oldPos(:,1,i);
    end
    
    %Update stress fiber positions realtive to cell centroid
    SFpos(3:4,Ibound) = SFpos(3:4,Ibound)-d*dir;
    [SFpos(1,Ibound),SFpos(2,Ibound)] = cart2pol(SFpos(3,Ibound),SFpos(4,Ibound));
    SFpos(1,:) = mod(SFpos(1,:),2*pi);
    aSF.pos(:,idx(i)+1:idx(i+1)) = SFpos;
    aSF.dir(:,i) = dir;
end
end

%**************************************************************************
%% Koff calculation for catch bonds
function k = kOff(f)
k=zeros(size(f));

kb = 1.3806e-23;
T = 310.15;
fx = 5.3800e-12;
xi = kb*T/fx;

%Rate constants
a = 3.309;
b = 0.0003942;
c = 58.19;

%This equation I made to fit the data from Kong 2009, JCB (Inversed to go from <t> to koff)
k(f>0) = (a*exp(-f(f>0)*xi./(kb*T))+(b*exp(f(f>0)*xi./(kb*T))+c*exp(-f(f>0)*xi./(kb*T))).^-1).^-1;
k(f<=0) = (a+(b+c).^-1).^-1;
end

%**************************************************************************
%% Force calculation based on substrate stiffness
function f=Force(Fs,t,pos,Kstiff,Ksoft,grad,v0)


%Shape of stiffness gradient
%  0      100 _____ 200          400 _____ 500
%            /     \                /     \
%           /       \              /       \
%          /         \            /         \
%_________/           \__________/           \________


wSoft = 170;%.................Width of soft substrate (um)
wStiff = 115;%................Width of stiff substrate (um)
w = grad;%.......................Width of stiffness gradient (um)

l1 = wSoft/2;
l2 = l1+wStiff;
p = wSoft+wStiff;


m = (Kstiff-Ksoft)/w;
b1 = Ksoft-m*(l1-w/2);
b2 = Kstiff+m*(l2-w/2);

x = pos;
Keff = zeros(size(x));

%Indexes to determine substrate location
Isoft = mod(abs(x),p)<l1-w/2 | mod(abs(x),p)>=l2+w/2;
Igrad1 = mod(abs(x),p)>=l1-w/2 & mod(abs(x),p)<l1+w/2;
Istiff = mod(abs(x),p)>=l1+w/2 & mod(abs(x),p)<l2-w/2;
Igrad2 = mod(abs(x),p)>=l2-w/2 & mod(abs(x),p)<l2+w/2;

%Set stiffness for position
Keff(Isoft) = Ksoft;
Keff(Istiff) = Kstiff;
Keff(Igrad1) = m.* mod(abs(x(Igrad1)),p)+b1;
Keff(Igrad2) = -m.* mod(abs(x(Igrad2)),p)+b2;

%Force calculation (2-spring model, Shwartz 2006, Biosystems)
f=Fs.*(1-exp(-(v0*(Keff.*t))./Fs));
end

%**************************************************************************
%% Speed and displacement calculation
function [avgVel,displ] = speed(rw,runTime)

    dim = size(rw,1);
    dsp = vecnorm(rw(1:dim,2:end)-rw(1:dim,1:end-1));
    dst = sum(dsp);%............................Gets distance cell traveled 
    avgVel = dst/runTime;
    displ = norm(rw(1:dim,end)-rw(1:dim,1));
end

%**************************************************************************
%% Fast binomial random number generator
function rnd = fbinornd(n,p)
    rnd = zeros(size(n));       
    I = find(n);     
    i = 1;
    if isscalar(p)
        while i<=numel(I)
            rnd(I(i)) = sum(rand(n(I(i)),1)<p);
            i=i+1;
        end
    else
        while i<=numel(I)
            rnd(I(i)) = sum(rand(n(I(i)),1)<p(I(i)));
            i=i+1;
        end
    end
end

%**************************************************************************
%% Fast normal random number generator
function rnd = fnormrnd(mu,sigma,varargin)
if ~numel(varargin)
    rnd = randn(1) .* sigma + mu;
else
    sizeOut=[varargin{:}];
    rnd = randn(sizeOut) .* sigma + mu;
end
end

%**************************************************************************
%% Fast poisson random number generator
function rnd = fpoissrnd(mu,varargin)

if ~numel(varargin)
    L = exp(-mu);
    p=1;
    k=-1;
    while p>L
        k=k+1;
        u = rand;
        p=p*u;
    end
    rnd = k;
else
    sizeOut=[varargin{:}];
    rnd = zeros(sizeOut);
    for i=1:numel(rnd)
        L = exp(-mu);
        p=1;
        k=-1;
        while p>L
            k=k+1;
            u = rand;
            p=p*u;
        end
        rnd(i) = k;
    end
end
rnd(rnd<0)=0;
end

%**************************************************************************
%% Display simulation time
function tclock=dToc(varargin)
tclock = [floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];

if numel(varargin) >0
    if isstring(varargin{1})
        title = varargin{1};
        fprintf('%s: %d:%0.2d:%0.5f \n',title,tclock);
    elseif varargin{1} == 0
        %Do nothing
    end
else
    fprintf('%d:%0.2d:%0.5f \n',tclock);
end
end

%**************************************************************************