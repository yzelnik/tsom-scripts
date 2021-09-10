function [stout,fext]=loadstfromp2p(p2pstr,Ps)

[xgrid,ygrid] = meshgrid((0.5:Ps.Nx)*Ps.Lx/Ps.Nx-Ps.Lx/2,(0.5:Ps.Ny)*Ps.Ly/Ps.Ny-Ps.Ly/2);

% try to get the number of points on the grid in x & y
dims = [length(unique(p2pstr.pdeo.grid.p(2,:))) length(unique(p2pstr.pdeo.grid.p(1,:)))];

xtmp=reshape(p2pstr.pdeo.grid.p(1,:),dims(1),dims(2));
ytmp=reshape(p2pstr.pdeo.grid.p(2,:),dims(1),dims(2));
tmpst = [];
for vind=1:Ps.VarNum
	data = reshape(p2pstr.u((1:p2pstr.np)+(vind-1)*p2pstr.np),dims(1),dims(2));
	tmpst(:,:,vind) = interp2(xtmp,ytmp,data,xgrid,ygrid)';
end;
stout=reshape(tmpst,Ps.Nx*Ps.Ny,Ps.VarNum);

% try to get the heterogenous f
par=p2pstr.u(p2pstr.nu+1:end); kf=0.15343097; 
tmpf=par(6)*(1+0.5*(1+cos((kf)*xtmp)));
fext = reshape(interp2(xtmp,ytmp,tmpf,xgrid,ygrid)',Ps.Nx*Ps.Ny,1);

end