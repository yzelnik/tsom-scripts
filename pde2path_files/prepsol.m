function p=prepsol(sizes,par,name,Vs,Ps)
p=zminit(sizes(1),sizes(2),sizes(3),par,name); 

tmpsol=repmat(reshape(Vs,Ps.Nx,Ps.Ny,3),1,2);
tmpsol=[tmpsol(:,end,:) tmpsol tmpsol(:,1,:)];
tmpsol=[tmpsol(end,:,:) ;tmpsol; tmpsol(1,:,:)];

tmpp=p; 
tmpp.nc.imax=20;
[xgrid,ygrid] = meshgrid((-0.5:Ps.Nx+1)*Ps.Lx/Ps.Nx-Ps.Lx/2,(-0.5:Ps.Ny*2+1)*Ps.Ly/Ps.Ny-Ps.Ly);

for ii=1:Ps.VarNum
    newintp = interp2(xgrid,ygrid,tmpsol(:,:,ii)',p.pdeo.grid.p(1,:),p.pdeo.grid.p(2,:));
    tmpp.u((1:tmpp.np)+(ii-1)*tmpp.np)=newintp;
end;
[nlout,i1,i2]=nloop(tmpp,tmpp.u);
disp(sprintf('%d steps, res: %.e.',i2,i1));
p.u=nlout; 

end