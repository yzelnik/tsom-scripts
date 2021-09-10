%% Scripts for the grazing management model, according to figure order
% 1) Bifurcation diagram and two phase-spaces (7 sections)
% 2) Snapshops of uniform vs. non-uniform grazing (2 sections)
% 3) Impact of temporal noise (2 sections)

%% setup bf data
% before running this section, run the AUTO script to create the 4 "b." and "s." files being read here (autoscript.txt)
filename = 'GrzMng_tur';
resratio = 4;
pntnum = 5000;

% Read uniform branches
unf = ReadAutoBif('b.GrzMng_unf',[1 2 -1 3 4],0,1);
bs1 = StabBfByLP(unf(1:pntnum,:),[1 find(unf(:,3)==1)]);
uv1 = StabBfByLP(unf(pntnum+1:pntnum*2,:),find(unf(pntnum+1:pntnum*2,3)==3));

% read periodic/Turing branch
tur = ReadAutoBif(['b.' filename],[1 2 -1 8]);
lx = resratio * tur(1,4);

% define model and other parameters
Ps=struct('LocFunc',@L_GrzMng,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'P',80,'eta',1.5,'kappa',1,'mu',2,'nu',5,'lambda',0.1,'gamma',4.,'rho',0.1,'Ds',[0.01 10],'VarNum',2,'Lx',lx,'Ly',1,'Nx',800,'Ny',1);
Es=struct('TsSize',0.1,'TimeDst',200,'JacMode',0,'VarInd',1,'LsaThresh',1e-4,'InitByOde',1,'StSmall',0.001,'SsThresh',1e-8,'St1Color',circshift(hsv(6),-2)*0.66,'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1],'BfStyle',['- ';'--';'x ']);

% add stability to periodic/Turing branch
turstb = AnalyzeAutoStates(filename,Ps,Es,@T_LSA,[],[1 2 8],{'P','','Lx'},'Es.MultPeriod',4);
% Try and run this -- if altpdist2/pdist2 is not available, skip this line and un-comment the next (manual)
tr1=RephaseBf(tur,turstb,Es,'Es.RephaseMode',1);
%tr1=[tur [zeros(347,1);ones(1520,1);zeros(728,1)]]; % manual fix

% plot this out, see what we got
plotbf({uv1,bs1,tr1},Es)
axis([70 140 -0.02 0.5])
%% split branches into parts to simplify later things

ratio = round(Ps.Lx/tur(1,4));

sts = ReadAutoStates(['s.' filename],Ps,Es,0,'Ps.Nx',Ps.Nx/ratio);

% splitting the uniform states
lps = sort([find(unf(:,3)==9) ; find(unf(:,3)==2)]);
clear bfunf; ind=1;
for ii=1:length(lps)-1
    tmp = unf(lps(ii)+1:lps(ii+1)-1,:);
    if(sum(tmp(:,2)<-1e-8)<10)
        bfunf{ind}=unf(lps(ii)+1:lps(ii+1)-1,:);
        ind=ind+1;
    end;
end;

% splitting the periodic states
tmp = tur(logical(tur(:,3)),1:3);
lps = [1 ; find(tmp(:,3)==5) ; size(tmp,1)];
clear stper; clear bfper;
for ii=1:length(lps)-1
    stper{ii}=sts(:,:,lps(ii):lps(ii+1));
    bfper{ii}=tmp(lps(ii):lps(ii+1),:);
end;

%% get specific points for the two phase-spaces

clear pnts;
maxnum = 20;
ds     = 0.02;
ppp    = [95 89];
porg   = 115;
orgind = 3;
tmax   = 200;
trjnum = 4;

% get the starting conditions
[~,loc]=min(abs(bfunf{orgind}(:,1)-porg));
orgst=NewtonLoop(repmat(bfunf{orgind}(loc,4:5),Ps.Nx,1),Ps,Es,'Ps.P',porg);

% For both weak and strong drought conditions
for pind= 1:2
  pval   = ppp(pind);
  ind    = 0;

  pnts=[];
  % calc uniform steady-state
  for ii=1:length(bfunf)
    if(pval<max(bfunf{ii}(:,1)) && pval>min(bfunf{ii}(:,1)))
        ind=ind+1;
        [~,loc]=min(abs(bfunf{ii}(:,1)-pval));

        pnts(:,:,ind)=NewtonLoop(repmat(bfunf{ii}(loc,4:5),Ps.Nx,1),Ps,Es,'Ps.P',pval);
    end;
  end;

  % calc periodic steady-state
  for ii=1:length(bfper)
    if(pval<max(bfper{ii}(:,1)) && pval>min(bfper{ii}(:,1)))
        ind=ind+1;
        [~,loc]=min(abs(bfper{ii}(:,1)-pval));

        [pnts(:,:,ind),loop,res]=NewtonLoop(repmat(stper{ii}(:,:,loc),round(ratio),1),Ps,Es,'Ps.P',pval);
    end;
  end;
  points{pind}=pnts;

end;
%% calc stable manifold (for weak drought scenario)
pind= 1;
pval= ppp(pind);
ind = 3;
sol = points{pind}(:,:,ind);
[stb,evs]=T_LSA(sol,Ps,Es,'Ps.P',pval);
%disp(stb)

[ee1,ee2]=eig(full(getjac(sol,Ps,Es)));
[val,ord]=sort(real(diag(ee2)));
actnum = min([maxnum sum(val>0)]);
steps = reshape(real(ee1(:,ord(end-actnum+1:end))),Ps.Nx,Ps.VarNum,actnum);

mannum=0;
for manind=1:actnum
	regnum = T_CountRegions(steps(:,:,manind),Ps,Es);
	if(regnum(1)==ratio) && (mannum<1)
        [frms,manifold] = runframes(sol+steps(:,:,manind)*ds,Ps,Es,'Ps.P',pval,'Es.Tstep',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
        mannum = mannum+1;
        if(manifold(1,2)>1e-3) && (manifold(1,3)>1e-3)
            [frms,manifold] = runframes(sol-steps(:,:,manind)*ds,Ps,Es,'Ps.P',pval,'Es.Tstep',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
        end;
    end;
end;

mans{pind}=manifold;

%% calc the unstable manifold (for strong drought scenario)
pind=2;
pval  = ppp(pind);
mults = [1 10 40]/1000;
reps  = 10;
mixind = size(pnts,3)-1;
twoends = [1 2];
twovals = mults([2 3]);
for ii=1:reps % binary search
    midval = mean(twovals);
    [frms,traj] = runframes(orgst+midval*pnts(:,:,mixind),Ps,Es,'Ps.P',pval,'Es.TsSize',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
    if(T_L2Norm(frms(:,:,end)-pnts(:,:,twoends(1)),Ps,Es)<T_L2Norm(frms(:,:,end)-pnts(:,:,twoends(2))))
        twovals(1)=midval;
    else
        twovals(2)=midval;
    end;
end;
% look how long the two ends don't diverge too far
[~,traj1] = runframes(orgst+twovals(1)*pnts(:,:,mixind),Ps,Es,'Ps.P',pval,'Es.Tstep',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
[~,traj2] = runframes(orgst+twovals(2)*pnts(:,:,mixind),Ps,Es,'Ps.P',pval,'Es.Tstep',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
manifold = traj(sqrt(sum((traj1(:,2:3)-traj2(:,2:3)).^2,2))<1e-3,:);
mans{pind}=manifold;

%% calc trajectories
% for both drought scenarios, calculate how the ecosystem responds to a drought, given its initial conditions
for pind=1:2
    pval   = ppp(pind);
    mixind = size(pnts,3)-1;

    % different intial conditions
    mults = [2 10 40]/1000;

    for ii=1:length(mults)
        [frms,traj] = runframes(orgst+mults(ii)*pnts(:,:,mixind),Ps,Es,'Ps.P',pval,'Es.Tstep',0.02,'Es.Frames',0:0.2:tmax,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
        tt{ii,pind}=traj;
    end;
end;


%% plot bifurcation diagram and two phase-spaces
specialpoints = [ppp(2) min(bfunf{2}(:,1)) ppp(1)];
nicestates=[frms(:,1,end),repmat([0.4 0],size(frms,1),1)];

% setup figure structure
clf;
axpos = [0.07 0.1 0.44 0.88;  0.62 0.1 0.36 0.39 ;  0.62 0.59 0.36 0.39];
bfcolors = [0 0 0; 0 0.5 0;0.2 0.8 0];
for ii = 1:size(axpos,1)
        ha(ii) = axes('Units','normalized','Position',axpos(ii,:));
end;
axes(ha(1))

% plot main bifurcation diagram
plotbf({bs1,uv1,tr1},Es,'Es.BfLineWidth',1.5,'Es.BfColor',bfcolors)
hold on;
plot(specialpoints(1)*[1 1],[-1 1],':k',specialpoints(2)*[1 1],[-1 1],':r',specialpoints(3)*[1 1],[-1 1],':k','lineWidth',1.5);
hold off;
axis([79 121 -0.02 0.44])

ylabel('biomass','fontSize',20);
xlabel('precipitation','fontSize',20);



epsi = [-0.002 0.00];
clrs = {repmat(bfcolors(3,:),3,1),[1 0 0; 1 0 0 ;bfcolors(3,:)]}; tinds = {2,1:2};
tinds = {2,[2 3]};

% plot each phase-space
for pind=1:2
    axes(ha(4-pind));

    % plot full and hollow circles in phase-spaces, according to their stability
    hold on;
    tmppnts = points{pind};
    for ii=1:size(tmppnts,3)
        sol = tmppnts(:,:,ii);
        [stb,evs]=T_LSA(sol,Ps,Es,'Ps.P',pval);
        tmp=T_MultiTest(sol,Ps,Es,'Es.TestFunc',{@T_CompFFT,@T_MinMax},'Es.CompFFT',[0 ratio]);
        if(stb)
            plot(tmp(:,1),tmp(:,2),'.k','markerSize',48);
        else
            plot(tmp(:,1),tmp(:,2),'ok','markerSize',14);
        end;
    end;

    % plot manifold
    tmpln = round(length(mans{pind})*0.9); % just to move lines slightly
    tmp = (tmpln:-1:1)'/tmpln*epsi(pind); tmp(length(mans{pind}))=0;  % just to move lines slightly
    plot(mans{pind}(:,2)+tmp,mans{pind}(:,3),'-b','lineWidth',2);

    % plot trajectories
    for ii=tinds{pind}
        tmpln = round(length(tt{ii,pind})*0.9);  % just to move lines slightly
        tmp = (tmpln:-1:1)'/tmpln*epsi(pind); tmp(length(tt{ii,pind}))=0;  % just to move lines slightly
        plot(tt{ii,pind}(:,2)-tmp,tt{ii,pind}(:,3),'-','color', clrs{pind}(ii,:),'lineWidth',1.5);
    end;
    hold off;
end;

ylabel('                           A_k (periodic mode amplitude)','fontSize',20);
xlabel('A_0 (uniform mode amplitude)','fontSize',20);

%% Simulate non-uniform grazing

Ps1d=struct('LocFunc',@L_GrzMng,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'P',80,'eta',1.5,'kappa',1,'mu',2,'nu',5,'lambda',0.1,'gamma',4.,'rho',0.1,'Ds',[0.01 10],'VarNum',2,'Lx',20,'Ly',1,'Nx',800,'Ny',1);
Es=struct('TsSize',0.05,'TimeDst',200,'JacMode',0,'VarInd',1,'InitByOde',1,'StSmall',0.01,'SsThresh',1e-6,'St1Color',circshift(hsv(6),-2)*0.66,'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1],'StAxis',[0 0.5]);

mshare  = 0.1; % how much of paramter M changes?
peropts = [16 8];
lenopts = round(Ps1d.Nx./peropts);

nugtime = 5;  % non-uniform grazing time
tottime = 50; % total simulatio time

% get to equilibrium at high P
st = run2ss(0.4,Ps1d,Es,'Es.OlDraw',0,'Es.TsSize',0.05,'Ps.P',115,'Es.InitFunc',@M_InitRndSt,'Es.StNoise',1e-4);

% scenario 1 -- keep uniform grazing throughout simulation
sts{1} = runframes(st,Ps1d,Es,'Es.OlDraw',0,'Ps.P',90,'Es.InitFunc',1e-4,'Es.Frames',0:tottime+1);

% run scenarios of non-uniform grazing
for ind=1:length(lenopts)
    % get mask for spatial structure of non-uniform grazing
    perlen = lenopts(ind);
    mask = repmat([ones(perlen,1) ;zeros(perlen,1)],ceil(Ps1d.Nx/(2*perlen)),1);
    mask = mask(1:Ps1d.Nx);

    % setup the non-uniform grazing for the first 5 years of the drought
    Es.DynVal={(1-mshare)*Ps1d.mu + 2*mask*mshare*Ps1d.mu ; mask*0+Ps1d.mu};
    Es.DynPrm='mu';
    Es.DynVal=Es.DynVal([ones(1,nugtime) ones(1,2+tottime-nugtime)*2]);

    sts{ind+1} = runframes(st(:,:,end),Ps1d,Es,'Es.OlDraw',0,'Ps.P',90,'Es.StNoise',1e-4,'Es.Frames',0:tottime+1);
end;

%% Now plot non-uniform grazing results
snapchs = [0 2 5 50]+1; % which snapshots to show (which years)?
% plot out snapshots of the different grazing scenarios
for ii=1:length(sts)
    for jj=1:length(snapchs)
        subplot(length(sts),length(snapchs),(ii-1)*length(snapchs)+jj);
        plotst(sts{ii},Ps1d,Es,snapchs(jj));
    end;
end;

%% Simulations under temporal noise
Ps1d=struct('LocFunc',@L_GrzMng,'SpaFunc',@S_RD,'IntegFunc',@I_PSRD,'P',80,'eta',1.5,'kappa',1,'mu',2,'nu',5,'lambda',0.1,'gamma',4.,'rho',0.1,'Ds',[0.01 10],'VarNum',2,'Lx',20,'Ly',1,'Nx',800,'Ny',1);
Es=struct('TsSize',0.1,'TimeDst',200,'JacMode',0,'VarInd',1,'InitByOde',1,'StSmall',0.01,'SsThresh',1e-6,'St1Color',circshift(hsv(6),-2)*0.66,'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1],'BfStyle',['- ';'--';'x ']);

% get patterned and uniform states
pat = run2ss([0;0.5],Ps1d,Es,'Es.OlDraw',1,'Ps.P',115,'Es.InitFunc',@M_InitPerSt,'Es.InitPrm',8);
unf = run2ss(0.3,Ps1d,Es,'Es.OlDraw',1,'Ps.P',115,'Es.InitFunc',@M_InitRndSt,'Es.StNoise',1e-4);

% These different tests will be used to figure out the final state
Es.TestList={@T_AvgVal,[1 1],@T_MinMax,[1 1],@T_CountRegions,[1 2]};

tottime = 100; % total simulation time
rndflag = 0; % 1 for Gamma distribution, 0 for Gaussian distribution
noiselvl= 0.02; %strength of stochasticity, one of: [0, 0.02, 0.1]

% NOTE: WITH THESE VALUES SIMULATION TAKES A LONG TIME.
% FOR QUICKER SIMULATIONS, REPLACE WITH COMMENTED-OUT LINES BELOW
randnum = 20;  % number of differnet randomization
pvals   = 115-(0:1:35); % different values of P
shares  = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]; % different shares of intial states

%randnum = 2;  % number of differnet randomization
%pvals   = 115-(0:5:35); % different values of P
%shares  = [0.001 0.01 0.1 1]; % different shares of intial states

% Go through different randomizations, initial conditions, and values of P
scores=[];
for rndind=1:randnum
    for ii=1:length(shares)
      disp([rndind ii])
      for jj=1:length(pvals)
        % starting conditions (mix of uniform and periodic state)
        start = shares(ii)*pat + (1-shares(ii))*unf;
        curp  = pvals(jj); % value of P
        rng(ii+jj*1e2+rndind*1e4); % randomization seed

        if(noiselvl>0) % Is there noise?
          if(rndflag) % Gamma distribution
            pdyn=gamrnd((1/noiselvl)^2,curp*noiselvl^2,tottime+1,1);
          else  % Gaussian distribution
              pdyn=max(0,randn(tottime+1,1)*noiselvl*curp+curp);
          end;
        else % no noise
          pdyn=repmat(curp,tottime+1,1);
        end;
          % run simulation
          tmpfrms = runframes(start,Ps1d,Es,'Es.DynVal',pdyn,'Es.DynPrm','P','Es.OlDraw',0,'Es.Frames',0:tottime);
          % save results of tests on final state from the simulation
          scores(ii,jj,:,rndind)=T_MultiTest(tmpfrms(:,:,end),Ps1d,Es);
      end;
    end;
end;

%% Plot out the results of the simulation grid
thresh=0.05;
finmat=squeeze(2*(scores(:,:,1,:,:)>thresh)-((scores(:,:,2,:,:)>thresh)));
imagesc(pvals,0.5:length(shares),flip(mean(finmat,3)));
