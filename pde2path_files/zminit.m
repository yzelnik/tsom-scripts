function p=zminit(lx,ly,h,par,name) 
p=[]; p=stanparam(p); pde=stanpdeo2D(lx,ly,h);
p.nc.tol=1e-6; p.nc.imax=50; p.file.smod=10; p=setfn(p,name); p.nc.ilam=9; 
p.nc.neq=3; p.plot.bpcmp=0; p.nc.lammin=0.5; p.nc.lammax=2; p.sw.sfem=-1; 
p.fuha.sG=@zmsG; p.fuha.sGjac=@zmsGjac; p.plot.pstyle=2; p.plot.axis='equal';
p.np=pde.grid.nPoints; x=pde.grid.p; p.pdeo=pde; 
screenlayout(p);
p.nu=p.np*p.nc.neq; p.sw.jac=1; 
uz=zeros(p.np,1); u=zeros(p.nu,1); icsw=1; % choose one 
switch icsw
  case 1; amp=0.0; w=0.2; x1=0; u(1:p.np)=uz; 
       u(p.np+1:2*p.np)=uz; u(2*p.np+1:3*p.np)=uz;
end
p.u=[u; par]; 
p=oosetfemops(p); 
p.nc.dsmin=0.00001; p.nc.dsmax=0.05;
p.sol.xi=0.1/p.nu; p.sw.spcalc=0; p.sw.bifcheck=1; p.nc.ntot=2000; 
plotsol(p,1,1,1); p.r=resi(p,p.u); fprintf('inires=%g\n',norm(p.r,Inf)); 
[p.u,p.r]=nloop(p,p.u); fprintf('res=%g\n',norm(p.r,Inf)); plotsol(p,1,1,1);  
p.nc.imax=10; 
p.sol.ds=0.001;
p.nc.neig=50; 
p.nc.dlammax=0.005; 
p.nc.nsteps=10;
p.sw.foldcheck=1;
p.sw.para=2;
%p.pm.mst=20;
p.file.smod=1; 

p.sw.spcalc=1;