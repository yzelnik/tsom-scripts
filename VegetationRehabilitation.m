%% Scripts for the vegetation rehabilitation model, according to figure order
% 1) Equilibrium and response to droughts (2 sections)
% 2) Impact of temporal noise (3 sections)
% 3) Simulations with parameters comparable to remote-sensing data (2 sections)

%% Initial setup (including perliminary simulations)

lamf=3.25; % baseline spatial scale
nop=5;  % number of periods
nx=128; dx=lamf*nop/nx;
ny=128; dy=lamf*nop/ny;

rng(1);
[xx,yy]=meshgrid((0.5:nx)*dx,(0.5:ny)*dy);
fnoise = 0.001;
kf = (nop*2*pi/(nx*dx));
Es=struct('TsSize',0.0002,'TimeDst',10,'SsThresh',1e-4,'NonNeg',1,'StSmall',0.01,'VarInd',1,'TsNum',500);
Ps=struct('LocFunc',@L_VegRhb,'SpaFunc',@S_VegRhb,'IntegFunc',@I_FDE,'kappa',5,'mu',3,'q',0.25,'nu',10,'lambda',0.1,'gamma',10,'alpha',100,'eta',0,'rho',0,'f',0.1,'P',240,'Ds',[0.02 2 2],'VarNum',3,'Lx',nx*dx,'Ly',ny*dy,'Nx',nx,'Ny',ny,'Bc',1);
%disp(sprintf('N %.2f, A %.2f, Q %.2f, G %.2f, P %.2f, dw %.2f, dh %.2f x%.3f',[Ps.nu/Ps.mu Ps.alpha/Ps.mu Ps.q/Ps.kappa Ps.gamma*Ps.kappa/Ps.mu  Ps.P*Ps.lambda/(Ps.mu*Ps.nu) Ps.Ds(2)/Ps.Ds(1) Ps.Ds(3)*Ps.nu/(Ps.Ds(1)*Ps.lambda) sqrt(Ps.Ds(1)/Ps.mu)]));

%
pval=270;
tsz=0.5*1e-4;
%Vs=rand(Ps.Nx,Ps.VarNum);
% modulation of f parameter
Ps.f = reshape(Ps.f*(1+0.5*(1+cos(kf*xx'))+fnoise*rand(size(xx))-fnoise/2),Ps.Nx*Ps.Ny,1);
Ps.Ly=Ps.Lx;
% Getting an similar initial state quickly, by using a reaction-diffusion model (without the H^2 term)
stripes0 = run2ss([Ps.f*1-0.1 ones(length(Ps.f),2)],Ps,Es,'Es.OlDraw',1,'Ps.IntegFunc',@I_FDSIMP,'Ps.SpaFunc',@S_RD,'Es.TsSize',0.05,'Ps.P',pval*1.1,'Es.TsNum',10,'Es.SsThresh',0.1);
% Now with the H^2 term
stripes1 = run2ss(stripes0,Ps,Es,'Es.OlDraw',1,'Es.TsSize',tsz,'Ps.P',pval);

%% Run differnet drought simulations

pps    = [240 220 240 220]; % values of P of weak and strong droughts
frgs   = [0 0 4 4]; % divide into how many parts
tmfrms = [0:5:40]; % save which snapshots?

for ii=1:length(pps)
  disp(ii);
  frgnum=frgs(ii);
  if(frgnum) % apply a rhombic state?
    % the x and y transects
    tmp1=repmat([ones(floor(ny/frgnum/2),1);zeros(ceil(ny/frgnum/2),1)],frgnum,1)-0.5;
    tmp2=repmat([ones(floor(nx/nop),1);zeros(ceil(nx/nop),1)],ceil((nop+1)/2),1)-0.5;
    offset=floor(ny/nop/2); % correct offset
    % make the mask to apply on the initial state to make it rohmbic-like
    mask=sign(tmp1.*(tmp2(offset+1:nx+offset)'))/2+1/2;
    % apply mask, reshape, and set into state structure
    newbiom = reshape(mask',nx*ny,1).*stripes1(:,1);
    newst=[newbiom stripes1(:,[2 3])];
  else
    newst=stripes1; % just use stripes
  end;
    %plotst(newst,Ps,Es); pause;
    % run simulation of 20 years, before drought
    frags1 = runsim(newst,Ps,Es,'Es.TimeDst',20,'Es.TsSize',tsz,'Ps.P',pval);
    % now apply drought, and save snapshots
    [frags2,tmpbf] = runframes(frags1,Ps,Es,'Es.OlDraw',2,'Es.TsSize',tsz,'Ps.P',pps(ii),'Es.Frames',tmfrms,'Es.TestFunc',@T_AvgVal);
    allfrms{ii}=frags2;
    %allbfs(:,:,ii)=tmpbf;
    %save nse_dim
end;

snapchs = [1 2 4 9]; % which snapshots to show (which years)?
% plot out snapshots of the different grazing scenarios
for ii=1:length(allfrms)
    for jj=1:length(snapchs)
        subplot(length(allfrms),length(snapchs),(ii-1)*length(snapchs)+jj);
        plotst(allfrms{ii},Ps,Es,snapchs(jj),'Es.StAxis',[0 0.5]);
    end;
end;

% For bif. diagram, run zmcmds.m file inside the pde2path_files subfolder (after installing pde2path: http://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/)
% For phase-spaces, run same parts as in the GrazingManagment.m file

%% Initial setup (including perliminary simulations) for noise sims

% Define x and y sizes and resolution
nx=100;
ny=100;
nop=4; % number of periods
Lx=13;
Ly=17;
dy=Ly/ny;
dx=Lx/nx;
kx=(2*pi*nop)/Lx;
ky=kx*0.75;

% define sinusodial pattern (to get rhombic state)
[xx,yy]=meshgrid((0.5:nx)*dx,(0.5:ny)*dy);
coscos=(cos(xx*kx/2).*cos(yy*ky))';
st2=[coscos(:)*0.1+0.1 ones(length(coscos(:)),2)];

Ps2d=struct('LocFunc',@L_VegRhb,'SpaFunc',@S_VegRhb,'IntegFunc',@I_FDE,'kappa',5,'mu',3,'q',0.25,'nu',10,'lambda',0.1,'gamma',10,'alpha',100,'eta',0,'rho',0,'f',0.1,'P',240,'Ds',[0.02 2 2],'VarNum',3,'Lx',Lx,'Ly',Ly,'Nx',nx,'Ny',ny,'Bc',1);
Es=struct('TsSize',0.00008,'TimeDst',200,'TimeMax',200,'JacMode',0,'VarInd',1,'InitByOde',1,'StSmall',0.01,'SsThresh',1e-5,'St1Color',circshift(hsv(6),-2)*0.66,'BfColor',[0 0 0 ; 1 0 0; 0 0.75 0; 0 0 1],'BfStyle',['- ';'--';'x ']);
fnoise=0;
Ps2d.f = reshape(Ps2d.f*(1+0.5*(1+cos(kx*xx'))+fnoise*rand(size(xx))-fnoise/2),Ps2d.Nx*Ps2d.Ny,1);
st1=[Ps2d.f*2-0.2 ones(length(Ps2d.f),2)];
plotst(st2,Ps2d,Es)
pval=260;

% short-cut using simple reaction-diffusion
altd = Ps2d.Ds; altd(3)=10;
stripes0 = run2ss(st1,Ps2d,Es,'Es.OlDraw',1,'Ps.IntegFunc',@I_FDSIMP,'Ps.SpaFunc',@S_RD,'Es.TsSize',0.05,'Ps.P',pval,'Ps.Ds',altd);
rhombic0 = run2ss(st2,Ps2d,Es,'Es.OlDraw',1,'Ps.IntegFunc',@I_FDSIMP,'Ps.SpaFunc',@S_RD,'Es.TsSize',0.05,'Ps.P',pval,'Ps.Ds',altd);
% and now the real thing
stripes = run2ss(stripes0,Ps2d,Es,'Es.OlDraw',1,'Ps.P',pval);
rhombic = run2ss(rhombic0,Ps2d,Es,'Es.OlDraw',1,'Ps.P',pval);

%% Now, simulations under drought with temporal noise

% These different tests will be used to figure out the final state
Es.TestList={@T_AvgVal,[1 1],@T_MinMax,[1 1],@T_CountRegions,[1 2]};

tottime = 100; % total simulation time
rndflag = 0; % 1 for Gamma distribution, 0 for Gaussian distribution
noiselvl= 0.02; %strength of stochasticity, one of: [0, 0.02, 0.1]

% NOTE: WITH THESE VALUES SIMULATION TAKES AN EXTREMELY LONG TIME.
% FOR SOMEWHAT-QUICKER SIMULATIONS, REPLACE WITH COMMENTED-OUT LINES BELOW
randnum = 20;  % number of differnet randomization
pvals   = 280-(0:2:60);% different values of P
shares  = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]; % different shares of intial states

%randnum = 2;  % number of differnet randomization
%pvals   = 280-(0:10:60); % different values of P
%shares  = [0.001 0.01 0.1 1]; % different shares of intial states

% Go through different randomizations, initial conditions, and values of P
scores=[];
for rndind=1:randnum
    for ii=1:length(shares)
      disp([rndind ii])
      for jj=1:length(pvals)
        % starting conditions (mix of uniform and periodic state)
        start = shares(ii)*rhombic + (1-shares(ii))*stripes;
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
          tmpfrms = runframes(start,Ps2d,Es,'Es.DynVal',pdyn,'Es.DynPrm','P','Es.OlDraw',1,'Es.Frames',0:tottime);
          % save results of tests on final state from the simulation
          scores(ii,jj,:,rndind)=T_MultiTest(tmpfrms(:,:,end),Ps2d,Es);
      end;
    end;
end;

%% Plot out the results of the simulation grid
thresh=0.01;
finmat=squeeze((scores(:,:,1,:,:)>thresh)+(scores(:,:,3,:,:)>1));
imagesc(pvals,0.5:length(shares),flip(mean(finmat,3)));

%% Simulations with different parameters (better matching data)

Lx=40; Ly=40; % system size
nx=100; ny=100; % resolution
Ps=struct('LocFunc',@L_VegRhb,'SpaFunc',@S_VegRhb,'IntegFunc',@I_FDE,'kappa',10,'mu',5,'q',1,'nu',7,'lambda',0.15,'gamma',5,'alpha',50,'eta',0,'rho',0,'f',0.15,'P',200,'Ds',[0.4 0.1 5],'VarNum',3,'Lx',Lx,'Ly',Ly,'Nx',nx,'Ny',ny,'Bc',1);
Es=struct('TsSize',0.0001,'TimeDst',10,'SsThresh',1e-5,'NonNeg',1,'StSmall',0.01,'VarInd',1,'TsNum',20,'TestFunc',@T_AvgVal);


rndind=19; % choose randomization of noise
kf=5*2*pi/Lx;
[xx,yy]=meshgrid((0.5:nx)*Lx/nx,(0.5:ny)*Ly/ny);

frmtms=[0 5 10 20 50 100 200];
rng(rndind)

% define the modulation of parameter f
fnoise=0.01;
fconst=0.1;
Ps.f = reshape(fconst*(1+0.5*(1+cos(kf*xx'))+fnoise*rand(size(xx))-fnoise/2),Ps.Nx*Ps.Ny,1);

% Run at high P (before drought)
tic;
out1 = runframes(repmat((Ps.f-0.1)*5,1,3),Ps,Es,'Ps.P',205,'Es.Frames',0:2:20,'Es.OlDraw',1);
toc;

% Add some noise to initial conditions
extnoise=0.1;
start=out1(:,:,end).*(1+rand(size(Ps.f))*extnoise-extnoise/2);
% Run at low P (after drought starts)
tic;
out2 = runframes(start,Ps,Es,'Ps.P',180,'Es.Frames',frmtms,'Es.OlDraw',1);
toc;
%% plot out results
fftsz=7; % number of frequencies in FFT
clf;
snapnum = size(out2,3);
for ii=1:snapnum
    cursnap=reshape(out2(:,1,ii),100,100);
    cursnap=(cursnap(11:90,11:90));
    % plot the state
    subplot(2,snapnum,ii)
    imagesc(cursnap');
    % now plot the FFT of the state
    subplot(2,snapnum,ii+snapnum)
    plotfft(cursnap,fftsz);
end;
