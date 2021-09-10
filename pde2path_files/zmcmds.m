%% start-up
close all; clear all; clear classes;
load rdm_sz64

par=[0.05; 10/3; 100/3; 0; 0; 0.1; 10000; 100; Ps.P; 50/3]; % Non-dimensional parameters
%    q    nu  al  eta   r    f  dh   dw   p  gam
lx=Ps.Lx/2; ly=Ps.Ly/2; h=1;

%% get bare-soil branch

bs1=prepsol([lx ly h],par,'bs1',ubsst,Ps);
%%
bs2=setfn(bs1,'bs2');
bs2.sol.ds=-0.001;

bs1=cont(bs1,5);
bs1=pmcont(bs1);
bs1=pmcont(bs1);

bs2=cont(bs2,5);
bs2=pmcont(bs2);

%% get stripes branch

str1=prepsol([lx ly h],par,'str1',strst,Ps);
str2=setfn(str1,'str2');
str2.sol.ds=-0.00001;
%
str1=cont(str1,5);
str1=pmcont(str1);
%str1=pmcont(str1);
%str1=cont(str1,5);
%str1=pmcont(str1);
%
str2=cont(str2,2); % used to be 5
str2=pmcont(str2);

str2=pmcont(str2);
%% get rhombic branch

rmb1=prepsol([lx ly h],par,'rmb1',rhost,Ps);
rmb1.pm.mst=20;

rmb2=setfn(rmb1,'rmb2');
rmb2.sol.ds=-0.001;

rmb1=cont(rmb1,20);
%rmb1=cont(rmb1,5);
rmb1=pmcont(rmb1);
%%

rmb2=cont(rmb2,5);
rmb2=pmcont(rmb2);
rmb2=pmcont(rmb2);
%% plot these branches
clf
plotbraf('bs1','pt90',3,0,'cl','k');
plotbraf('bs2','pt39',3,0,'cl','k');
plotbraf('rmb1','pt150',3,0,'cl','r');
plotbraf('rmb2','pt180',3,0,'cl','r');
plotbraf('str1','pt90',3,0,'cl','b');
plotbraf('str2','pt150',3,0,'cl','b');

xlabel('p'); ylabel('||b||');
axis([0.65 1.3 0 6]);

hold on; plot([0.8 0.8],[0 10],'k:'); hold off;
%ha=axes('Position',[.11 .58 .3 .3]); box on;
%plotsolf('rmb1','pt50',0,1,2); title('rhombic');
%ha=axes('Position',[.64 .38 .3 .3]); box on;
%plotsolf('str1','pt50',0,1,2); title('stripe');
%% load pt from rmb2, meshrefine, then pmcont; -works!
p=loadp('rmb2','pt5'); plotsol(p);
p.fuha.e2rs=@sgeefu; p=oomeshada(p,'ngen',3,'sig',0.2); p.np
%% mesh-adaption (makes, e.g., hex-cont more reliable by making a 'hex-mesh')
%p=loadp('rmb2','pt10'); plotsol(p);
p=prepsol([lx ly h],par,'rmb',rhost,Ps);

p.fuha.e2rs=@sgeefu; p.nc.imax=20;
p=cont(p,1);
p=oomeshada(p,'ngen',4,'sig',0.1); p.np
p=cont(p,1);
p=oomeshada(p,'ngen',4,'sig',0.1); p.np
[p.mesh.bp,p.mesh.bt,p.mesh.be]=getpte(p);
%% now continue
p=setfn(p,'rmb2');
p.fuha.blss=@mbel; p.hopf.ilss=1; p.nc.mbw=2; p.sw.verb=2; p.sw.spcalc=1;
p.nc.amod=5; p.nc.ngen=4; p.nc.dsmax=0.002; p.sw.bifcheck=0; p.sw.foldcheck=0;
p.sol.ds=-1e-2; p.nc.tol=1e-8; p.nc.dlammax=0.01;
p=cont(p,100);
%p=pmcont(p);
%p=pmcont(p);
%% branch plot
figure(3); clf; cmp=0;% cmp=11;
plotbraf('rmb2','pt40',3,cmp,'cl','k','labi',10);
plotbraf('rmb2b','pt70',3,cmp,'cl','r','labi',20);
%% soln plot
for i=10:10:40; plotsol('rmb2',['pt' mat2str(i)]); axis image; pause; end
%%
for i=10:20:70; plotsol('rmb2b',['pt' mat2str(i)]); axis image; pause; end
%% other direction, still problematic beyond 0.912
p=loadp('rmb2','pt5'); plotsol(p); p=setfn(p,'rmb2b');
p.fuha.e2rs=@sgeefu; p.nc.tol=1e-8;
p=oomeshada(p,'ngen',2,'sig',0.3); p.np
[p.mesh.bp,p.mesh.bt,p.mesh.be]=getpte(p);
%%
p.fuha.blss=@mbel; p.hopf.ilss=1; p.nc.mbw=2; p.sw.verb=2; p.sw.spcalc=1;
p.nc.amod=5; p.nc.ngen=2; p.nc.dsmax=0.5; p.sw.bifcheck=0; p.sw.foldcheck=0;
p.sol.ds=-p.sol.ds; p.nc.dsmax=0.05; p=cont(p,30);
%%
figure(3); clf; cmp=0;% cmp=11;
plotbraf('rmb1','pt150',3,cmp,'cl','r','labi',20);
plotbraf('rmb2','pt180',3,cmp,'cl','r','labi',20);
plotbraf('str1','pt90',3,cmp,'cl','b','labi',20);
plotbraf('str2','pt170',3,cmp,'cl','b','labi',20);
